use crate::config::QcConfig;
use crate::io::io::open_fastq;
use crate::progress::progress::{ProgressTracker, QC_STAGE};
use anyhow::anyhow;
use flate2::Compression;
use flate2::write::GzEncoder;
use seq_io::fastq::Record;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::LazyLock;
use std::sync::atomic::{AtomicUsize, Ordering};

const TOTAL_IDX: usize = 0;
const KEPT_IDX: usize = 1;
const DROPPED_IDX: usize = 2;

/// FASTQ Phred ASCII offset (Sanger / Illumina 1.8+ / ONT).
const PHRED_OFFSET: u8 = 33;

/// `10^(-q/10)` for every possible raw Phred score, so scoring a read costs a table
/// lookup per base rather than a `powf`.
static ERROR_PROB: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0f64; 256];
    for (q, p) in table.iter_mut().enumerate() {
        *p = 10f64.powf(-(q as f64) / 10.0);
    }
    table
});

/// The drop tallies below are process-wide, so any test that asserts on them has to
/// hold this lock for the duration of its run.
#[cfg(test)]
pub(crate) static QC_COUNT_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Take [`QC_COUNT_TEST_LOCK`] and start from zeroed tallies, ignoring poisoning from
/// an unrelated failing test.
#[cfg(test)]
pub(crate) fn lock_qc_counts_for_test() -> std::sync::MutexGuard<'static, ()> {
    let guard = QC_COUNT_TEST_LOCK
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    reset_qc_counts();
    guard
}

static DROPPED_LENGTH: AtomicUsize = AtomicUsize::new(0);
static DROPPED_QUALITY: AtomicUsize = AtomicUsize::new(0);
static MISSING_QUALITY: AtomicUsize = AtomicUsize::new(0);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct QcCounts {
    pub dropped_length: usize,
    pub dropped_quality: usize,
    /// Reads seen while a quality threshold was active but that carried no quality string.
    pub missing_quality: usize,
}

pub fn qc_counts() -> QcCounts {
    QcCounts {
        dropped_length: DROPPED_LENGTH.load(Ordering::Relaxed),
        dropped_quality: DROPPED_QUALITY.load(Ordering::Relaxed),
        missing_quality: MISSING_QUALITY.load(Ordering::Relaxed),
    }
}

pub fn reset_qc_counts() {
    DROPPED_LENGTH.store(0, Ordering::Relaxed);
    DROPPED_QUALITY.store(0, Ordering::Relaxed);
    MISSING_QUALITY.store(0, Ordering::Relaxed);
}

/// Mean-error-probability read quality: `Q = -10 * log10(mean(10^(-q_i/10)))`.
///
/// This is the convention used by Dorado, NanoFilt, chopper and filtlong. It is
/// deliberately *not* the arithmetic mean of the Phred scores, which overestimates
/// quality whenever a read carries a few very bad bases.
///
/// Expects raw Phred values (no ASCII offset). Returns `None` for an empty input.
pub fn mean_error_quality<I: IntoIterator<Item = u8>>(phred: I) -> Option<f64> {
    let table = &*ERROR_PROB;
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for q in phred {
        sum += table[q as usize];
        n += 1;
    }
    if n == 0 {
        return None;
    }
    let mean_p = sum / n as f64;
    // Only underflows for absurdly high Phred values; treat those as perfect.
    if mean_p <= 0.0 {
        return Some(f64::INFINITY);
    }
    Some(-10.0 * mean_p.log10())
}

/// [`mean_error_quality`] for a FASTQ quality line (ASCII, +33 offset).
pub fn mean_error_quality_ascii(qual: &[u8]) -> Option<f64> {
    mean_error_quality(qual.iter().map(|&q| q.saturating_sub(PHRED_OFFSET)))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QcVerdict {
    Pass,
    FailLength,
    FailQuality,
}

/// Length / quality thresholds applied to a raw read. Used both standalone
/// (`sarracenia qc`) and as a pre-alignment gate inside `annotate` / `kit`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct QcFilter {
    pub min_length: Option<usize>,
    pub max_length: Option<usize>,
    pub min_quality: Option<f64>,
}

impl QcFilter {
    pub fn new(
        min_length: Option<usize>,
        max_length: Option<usize>,
        min_quality: Option<f64>,
    ) -> Self {
        Self {
            min_length,
            max_length,
            min_quality,
        }
    }

    /// True when at least one threshold is set; callers skip all QC work otherwise.
    pub fn is_active(&self) -> bool {
        self.min_length.is_some() || self.max_length.is_some() || self.min_quality.is_some()
    }

    /// True when a threshold requires actually scoring the quality string.
    pub fn needs_quality(&self) -> bool {
        self.min_quality.is_some()
    }

    /// Catch thresholds that can never be satisfied up front, rather than letting the
    /// user discover it from an empty output file. Quotes the pre-alignment flag names;
    /// use [`QcFilter::validate_flags`] for the post-trim gate.
    pub fn validate(&self) -> anyhow::Result<()> {
        self.validate_flags("")
    }

    /// [`QcFilter::validate`], but `flag_infix` is spliced into the flag names quoted in
    /// the error so it names the flag the user actually typed: `""` gives
    /// `--min-length`, `"trimmed-"` gives `--min-trimmed-length`.
    pub fn validate_flags(&self, flag_infix: &str) -> anyhow::Result<()> {
        if let (Some(min), Some(max)) = (self.min_length, self.max_length)
            && min > max
        {
            return Err(anyhow!(
                "--min-{flag_infix}length ({min}) is greater than --max-{flag_infix}length ({max}), no read can pass"
            ));
        }
        if let Some(q) = self.min_quality
            && (q.is_nan() || q < 0.0)
        {
            return Err(anyhow!(
                "--min-{flag_infix}quality must be a non-negative number, got {q}"
            ));
        }
        Ok(())
    }

    /// `quality` is the pre-computed mean-error quality, or `None` when the record
    /// carries no quality string. A read with no quality *passes* the quality check
    /// and is tallied in [`QcCounts::missing_quality`] rather than dropped silently.
    pub fn evaluate(&self, len: usize, quality: Option<f64>) -> QcVerdict {
        if let Some(min) = self.min_length
            && len < min
        {
            DROPPED_LENGTH.fetch_add(1, Ordering::Relaxed);
            return QcVerdict::FailLength;
        }
        if let Some(max) = self.max_length
            && len > max
        {
            DROPPED_LENGTH.fetch_add(1, Ordering::Relaxed);
            return QcVerdict::FailLength;
        }
        if let Some(min) = self.min_quality {
            match quality {
                Some(q) if q < min => {
                    DROPPED_QUALITY.fetch_add(1, Ordering::Relaxed);
                    return QcVerdict::FailQuality;
                }
                None => {
                    MISSING_QUALITY.fetch_add(1, Ordering::Relaxed);
                }
                _ => {}
            }
        }
        QcVerdict::Pass
    }

    pub fn passes(&self, len: usize, quality: Option<f64>) -> bool {
        self.evaluate(len, quality) == QcVerdict::Pass
    }

    /// Score a FASTQ quality line, but only if a threshold actually needs it.
    pub fn quality_of_ascii(&self, qual: &[u8]) -> Option<f64> {
        if self.needs_quality() {
            mean_error_quality_ascii(qual)
        } else {
            None
        }
    }

    /// Score raw Phred scores (BAM), but only if a threshold actually needs it.
    pub fn quality_of_phred<I: IntoIterator<Item = u8>>(&self, phred: I) -> Option<f64> {
        if self.needs_quality() {
            mean_error_quality(phred)
        } else {
            None
        }
    }

    /// Score and evaluate a read in one step; `qual` is a FASTQ quality line (ASCII, +33).
    pub fn passes_ascii(&self, seq_len: usize, qual: &[u8]) -> bool {
        self.passes(seq_len, self.quality_of_ascii(qual))
    }

    /// Score and evaluate a read in one step; `qual` is raw Phred (BAM).
    pub fn passes_phred(&self, seq_len: usize, qual: &[u8]) -> bool {
        self.passes(seq_len, self.quality_of_phred(qual.iter().copied()))
    }

    /// One-line human summary of the active thresholds, for the pipeline banner.
    pub fn describe(&self) -> String {
        let mut parts = Vec::new();
        if let Some(v) = self.min_length {
            parts.push(format!("min length {v}"));
        }
        if let Some(v) = self.max_length {
            parts.push(format!("max length {v}"));
        }
        if let Some(v) = self.min_quality {
            parts.push(format!("min quality Q{v}"));
        }
        parts.join(", ")
    }
}

/// Report reads dropped by a QC gate, plus reads the quality filter could not be
/// applied to. No-op when nothing was filtered.
pub fn print_qc_summary(qc: &QcFilter) {
    print_qc_summary_labeled(qc, "QC");
}

/// [`print_qc_summary`] with a caller-supplied label, so the pre-alignment gate and
/// the post-trim gate are distinguishable in a single `kit` run.
pub fn print_qc_summary_labeled(qc: &QcFilter, label: &str) {
    let counts = qc_counts();
    let dropped = counts.dropped_length + counts.dropped_quality;
    if dropped > 0 {
        println!(
            "{label} dropped {dropped} read(s) ({} by length, {} by quality)",
            counts.dropped_length, counts.dropped_quality
        );
    }
    if qc.needs_quality() && counts.missing_quality > 0 {
        println!(
            "Warning: {} read(s) carried no quality scores; the quality filter was not applied to them.",
            counts.missing_quality
        );
    }
}

/// FASTQ output sink that keeps the gzip encoder concrete so the trailer can be
/// written (and its errors surfaced) via an explicit `finish`.
enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl FastqWriter {
    fn create(path: &str, gzip: bool) -> anyhow::Result<Self> {
        let file = File::create(path)
            .map_err(|e| anyhow!("Failed to create output file '{path}': {e}"))?;
        let buffered = BufWriter::new(file);
        Ok(if gzip {
            FastqWriter::Gzip(GzEncoder::new(buffered, Compression::default()))
        } else {
            FastqWriter::Plain(buffered)
        })
    }

    fn write_record(&mut self, head: &[u8], seq: &[u8], qual: &[u8]) -> std::io::Result<()> {
        let w: &mut dyn Write = match self {
            FastqWriter::Plain(w) => w,
            FastqWriter::Gzip(w) => w,
        };
        w.write_all(b"@")?;
        w.write_all(head)?;
        w.write_all(b"\n")?;
        w.write_all(seq)?;
        w.write_all(b"\n+\n")?;
        w.write_all(qual)?;
        w.write_all(b"\n")
    }

    fn finish(self) -> anyhow::Result<()> {
        match self {
            FastqWriter::Plain(mut w) => w.flush()?,
            FastqWriter::Gzip(w) => {
                w.finish()?.flush()?;
            }
        }
        Ok(())
    }
}

fn wants_gzip(path: &str, flag: bool) -> bool {
    flag || path.to_ascii_lowercase().ends_with(".gz")
}

fn new_progress(output: &str, verbose: bool) -> ProgressTracker {
    if verbose {
        let log_dir = Path::new(output).parent().unwrap_or_else(|| Path::new("."));
        ProgressTracker::new_with_logging(&QC_STAGE, "qc", log_dir)
    } else {
        ProgressTracker::new(&QC_STAGE)
    }
}

/// Standalone length/quality filtering: reads in, surviving reads out, in the same
/// format as the input. Does no alignment.
pub fn run_qc(
    input: &str,
    output: &str,
    dropped: Option<&str>,
    config: &QcConfig,
) -> anyhow::Result<()> {
    config.filter.validate()?;
    if !config.filter.is_active() {
        return Err(anyhow!(
            "No thresholds given; set at least one of --min-length, --max-length, --min-quality"
        ));
    }
    reset_qc_counts();

    if input.to_ascii_lowercase().ends_with(".bam") {
        run_qc_bam(input, output, dropped, config)
    } else {
        run_qc_fastq(input, output, dropped, config)
    }
}

fn run_qc_fastq(
    input: &str,
    output: &str,
    dropped: Option<&str>,
    config: &QcConfig,
) -> anyhow::Result<()> {
    let qc = &config.filter;
    let mut reader = open_fastq(input);

    let mut kept_writer = FastqWriter::create(output, wants_gzip(output, config.gzip))?;
    let mut dropped_writer = match dropped {
        Some(path) => Some(FastqWriter::create(path, wants_gzip(path, config.gzip))?),
        None => None,
    };

    let progress = new_progress(output, config.verbose);

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| anyhow!("Input FASTQ parsing failed: {e}"))?;
        let seq = record.seq();
        let qual = record.qual();
        let quality = qc.quality_of_ascii(qual);

        progress.inc(TOTAL_IDX);
        if qc.passes(seq.len(), quality) {
            progress.inc(KEPT_IDX);
            kept_writer
                .write_record(record.head(), seq, qual)
                .map_err(|e| anyhow!("Failed to write to '{output}': {e}"))?;
        } else {
            progress.inc(DROPPED_IDX);
            if let Some(w) = dropped_writer.as_mut() {
                w.write_record(record.head(), seq, qual)
                    .map_err(|e| anyhow!("Failed to write dropped read: {e}"))?;
            }
        }

        progress.refresh();
    }

    kept_writer.finish()?;
    if let Some(w) = dropped_writer {
        w.finish()?;
    }
    progress.finish();
    print_qc_summary(qc);
    Ok(())
}

fn run_qc_bam(
    input: &str,
    output: &str,
    dropped: Option<&str>,
    config: &QcConfig,
) -> anyhow::Result<()> {
    use noodles::sam::alignment::io::Write as SamWrite;

    let qc = &config.filter;
    let mut reader = File::open(input)
        .map(noodles::bam::io::Reader::new)
        .map_err(|e| anyhow!("Failed to open BAM file '{input}': {e}"))?;
    // Reads are passed through untouched, so the input header stays valid.
    let header = reader.read_header()?;

    let mut kept_writer = noodles::bam::io::Writer::new(
        File::create(output).map_err(|e| anyhow!("Failed to create output file '{output}': {e}"))?,
    );
    kept_writer.write_header(&header)?;

    let mut dropped_writer = match dropped {
        Some(path) => {
            let mut w = noodles::bam::io::Writer::new(
                File::create(path)
                    .map_err(|e| anyhow!("Failed to create dropped file '{path}': {e}"))?,
            );
            w.write_header(&header)?;
            Some(w)
        }
        None => None,
    };

    let progress = new_progress(output, config.verbose);

    for result in reader.records() {
        let record = result.map_err(|e| anyhow!("Failed to read BAM record: {e}"))?;
        let len = record.sequence().len();
        let quality = qc.quality_of_phred(record.quality_scores().iter());

        progress.inc(TOTAL_IDX);
        if qc.passes(len, quality) {
            progress.inc(KEPT_IDX);
            kept_writer.write_alignment_record(&header, &record)?;
        } else {
            progress.inc(DROPPED_IDX);
            if let Some(w) = dropped_writer.as_mut() {
                w.write_alignment_record(&header, &record)?;
            }
        }

        progress.refresh();
    }

    progress.finish();
    print_qc_summary(qc);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn approx(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-6, "expected {b}, got {a}");
    }

    #[test]
    fn test_mean_error_quality_uniform() {
        // All Q40 bases -> mean error prob 1e-4 -> Q40 back out.
        approx(mean_error_quality_ascii(b"IIIIIIII").unwrap(), 40.0);
    }

    #[test]
    fn test_mean_error_quality_is_not_arithmetic_mean() {
        // Four Q40 bases and one Q0 base. The arithmetic mean would be 32.0; the
        // error-probability mean is dominated by the single bad base.
        let q = mean_error_quality_ascii(b"IIII!").unwrap();
        let expected_mean_p: f64 = ((4.0 * 1e-4) + 1.0) / 5.0;
        approx(q, -10.0 * expected_mean_p.log10());
        assert!(q < 7.0, "expected the bad base to dominate, got {q}");
    }

    #[test]
    fn test_mean_error_quality_empty_is_none() {
        assert!(mean_error_quality_ascii(b"").is_none());
    }

    #[test]
    fn test_mean_error_quality_raw_phred_matches_ascii() {
        let ascii = mean_error_quality_ascii(b"5?I").unwrap();
        let raw = mean_error_quality([b'5' - 33, b'?' - 33, b'I' - 33]).unwrap();
        approx(ascii, raw);
    }

    #[test]
    fn test_length_bounds() {
        let qc = QcFilter::new(Some(100), Some(200), None);
        assert_eq!(qc.evaluate(99, None), QcVerdict::FailLength);
        assert_eq!(qc.evaluate(100, None), QcVerdict::Pass);
        assert_eq!(qc.evaluate(200, None), QcVerdict::Pass);
        assert_eq!(qc.evaluate(201, None), QcVerdict::FailLength);
    }

    #[test]
    fn test_quality_threshold() {
        let qc = QcFilter::new(None, None, Some(10.0));
        assert_eq!(qc.evaluate(500, Some(9.99)), QcVerdict::FailQuality);
        assert_eq!(qc.evaluate(500, Some(10.0)), QcVerdict::Pass);
    }

    #[test]
    fn test_missing_quality_passes_and_is_counted() {
        let _guard = lock_qc_counts_for_test();
        let qc = QcFilter::new(None, None, Some(10.0));
        assert_eq!(qc.evaluate(500, None), QcVerdict::Pass);
        assert_eq!(qc_counts().missing_quality, 1);
    }

    #[test]
    fn test_inactive_filter_does_no_work() {
        let qc = QcFilter::default();
        assert!(!qc.is_active());
        assert!(!qc.needs_quality());
        assert!(qc.quality_of_ascii(b"IIII").is_none());
        assert_eq!(qc.evaluate(0, None), QcVerdict::Pass);
    }

    #[test]
    fn test_validate_rejects_impossible_bounds() {
        assert!(QcFilter::new(Some(500), Some(100), None).validate().is_err());
        assert!(QcFilter::new(Some(100), Some(500), None).validate().is_ok());
        assert!(QcFilter::new(None, None, Some(-1.0)).validate().is_err());
    }

    #[test]
    fn test_run_qc_fastq_end_to_end() {
        let _guard = lock_qc_counts_for_test();
        // read1: 8bp Q40 -> kept. read2: 4bp -> too short. read3: 8bp but Q0 -> too low.
        let content = b"@read1 desc\nACGTACGT\n+\nIIIIIIII\n\
                        @read2\nACGT\n+\nIIII\n\
                        @read3\nTTTTAAAA\n+\n!!!!!!!!\n";
        let mut input = NamedTempFile::with_suffix(".fastq").unwrap();
        input.write_all(content).unwrap();
        input.flush().unwrap();

        let out = NamedTempFile::with_suffix(".fastq").unwrap();
        let rejected = NamedTempFile::with_suffix(".fastq").unwrap();

        let config = QcConfig {
            filter: QcFilter::new(Some(8), None, Some(10.0)),
            gzip: false,
            verbose: false,
        };
        run_qc(
            input.path().to_str().unwrap(),
            out.path().to_str().unwrap(),
            Some(rejected.path().to_str().unwrap()),
            &config,
        )
        .unwrap();

        let kept = std::fs::read_to_string(out.path()).unwrap();
        // The full header (id + description) must survive the round trip.
        assert!(kept.contains("@read1 desc"), "got: {kept}");
        assert!(!kept.contains("@read2"));
        assert!(!kept.contains("@read3"));

        let dropped = std::fs::read_to_string(rejected.path()).unwrap();
        assert!(dropped.contains("@read2"));
        assert!(dropped.contains("@read3"));

        let counts = qc_counts();
        assert_eq!(counts.dropped_length, 1);
        assert_eq!(counts.dropped_quality, 1);
    }

    #[test]
    fn test_run_qc_rejects_empty_filter() {
        let input = NamedTempFile::with_suffix(".fastq").unwrap();
        let out = NamedTempFile::with_suffix(".fastq").unwrap();
        let config = QcConfig {
            filter: QcFilter::default(),
            gzip: false,
            verbose: false,
        };
        assert!(
            run_qc(
                input.path().to_str().unwrap(),
                out.path().to_str().unwrap(),
                None,
                &config,
            )
            .is_err()
        );
    }
}
