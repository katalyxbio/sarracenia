use std::fs::File;
use std::io::{BufWriter, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Suppresses every informational and progress line (`--quiet`). Errors are still
/// reported.
static QUIET: AtomicBool = AtomicBool::new(false);

pub fn set_quiet(quiet: bool) {
    QUIET.store(quiet, Ordering::SeqCst);
}

pub fn is_quiet() -> bool {
    QUIET.load(Ordering::SeqCst)
}

/// Print one informational line to stdout, unless `--quiet` is set.
pub fn info(line: impl std::fmt::Display) {
    if !is_quiet() {
        println!("{line}");
    }
}

/// How often the progress line is redrawn. On a terminal it is rewritten in place, so
/// it can be fast; piped into a file or CI log every redraw is a new line, so it is
/// slow enough that a long run leaves a handful of lines rather than thousands.
const TTY_REDRAW_MS: u64 = 100;
const PIPE_REDRAW_MS: u64 = 10_000;

/// Column the counters start at, so stage names line up across stages.
const NAME_WIDTH: usize = 12;

/// One stage's progress line. `metrics[0]` is always the running total (rendered as
/// `<total> <unit>`); the rest are rendered as `| <name> <count>`.
pub(crate) struct StageSpec {
    pub name: &'static str,
    pub unit: &'static str,
    pub metrics: &'static [&'static str],
}

pub(crate) const ANNOTATE_STAGE: StageSpec = StageSpec {
    name: "Annotating",
    unit: "reads",
    metrics: &["total", "kept", "dropped"],
};

pub(crate) const FILTER_STAGE: StageSpec = StageSpec {
    name: "Filtering",
    unit: "reads",
    metrics: &["total", "kept", "dropped"],
};

pub(crate) const QC_STAGE: StageSpec = StageSpec {
    name: "QC",
    unit: "reads",
    metrics: &["total", "kept", "dropped"],
};

pub(crate) const TRIM_STAGE: StageSpec = StageSpec {
    name: "Trimming",
    unit: "reads",
    metrics: &["total", "kept", "kept split", "failed"],
};

/// `"kept split"` -> `"Kept split:"`. Keeps the `--verbose` log file's metric column
/// in the format it has always had.
fn log_label(metric: &str) -> String {
    let mut chars = metric.chars();
    match chars.next() {
        Some(first) => format!("{}{}:", first.to_ascii_uppercase(), chars.as_str()),
        None => ":".to_string(),
    }
}

/// Draws the single progress line for a stage, either rewriting it in place (terminal)
/// or appending it (anything else).
struct Renderer {
    tty: bool,
    redraw_ms: u64,
    start: Instant,
    last_draw_ms: AtomicU64,
    /// Width of the line currently sitting on the terminal, so it can be blanked out
    /// before a shorter line replaces it.
    on_screen: Mutex<usize>,
}

impl Renderer {
    fn new() -> Self {
        let tty = std::io::stderr().is_terminal();
        Self {
            tty,
            redraw_ms: if tty { TTY_REDRAW_MS } else { PIPE_REDRAW_MS },
            start: Instant::now(),
            last_draw_ms: AtomicU64::new(0),
            on_screen: Mutex::new(0),
        }
    }

    /// True at most once per redraw interval, for exactly one caller — workers call
    /// `refresh` concurrently and must not all draw the same line.
    fn due(&self) -> bool {
        let now = self.start.elapsed().as_millis() as u64;
        let last = self.last_draw_ms.load(Ordering::Relaxed);
        if now.saturating_sub(last) < self.redraw_ms {
            return false;
        }
        self.last_draw_ms
            .compare_exchange(last, now, Ordering::Relaxed, Ordering::Relaxed)
            .is_ok()
    }

    fn draw(&self, line: &str, last: bool) {
        let mut on_screen = self
            .on_screen
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let mut err = std::io::stderr().lock();

        if self.tty {
            let blank = on_screen.saturating_sub(line.len());
            let _ = write!(err, "\r{line}{:blank$}", "");
            if last {
                let _ = writeln!(err);
                *on_screen = 0;
            } else {
                *on_screen = line.len();
            }
        } else {
            let _ = writeln!(err, "{line}");
            *on_screen = 0;
        }
        let _ = err.flush();
    }

    /// Blank the in-place line so a message can be printed without landing on top of it.
    fn erase(&self) {
        let mut on_screen = self
            .on_screen
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if self.tty && *on_screen > 0 {
            let mut err = std::io::stderr().lock();
            let _ = write!(err, "\r{:width$}\r", "", width = *on_screen);
            let _ = err.flush();
        }
        *on_screen = 0;
    }
}

pub(crate) struct ProgressTracker {
    spec: &'static StageSpec,
    counts: Vec<AtomicUsize>,
    error_msg: Mutex<Option<String>>,
    log: Option<ProgressLog>,
    /// `None` under `--quiet`, which is what makes every draw call a no-op.
    renderer: Option<Renderer>,
}

struct ProgressLog {
    path: PathBuf,
    step: String,
}

impl ProgressLog {
    fn write(&self, counts: &[AtomicUsize], spec: &StageSpec) -> Result<(), String> {
        let file = File::create(&self.path)
            .map_err(|e| format!("Failed to create log file '{}': {e}", self.path.display()))?;

        let mut w = BufWriter::new(file);
        writeln!(w, "step\tmetric\tcount")
            .map_err(|_| format!("Failed to write progress log '{}'", self.path.display()))?;

        for (count, metric) in counts.iter().zip(spec.metrics.iter()) {
            writeln!(
                w,
                "{}\t{}\t{}",
                self.step,
                log_label(metric),
                count.load(Ordering::Relaxed)
            )
            .map_err(|_| format!("Failed to write progress log '{}'", self.path.display()))?;
        }

        w.flush()
            .map_err(|_| format!("Failed to flush progress log '{}'", self.path.display()))
    }
}

impl ProgressTracker {
    pub(crate) fn new(spec: &'static StageSpec) -> Self {
        Self::new_inner(spec, None)
    }

    pub(crate) fn new_with_logging(
        spec: &'static StageSpec,
        step: impl Into<String>,
        log_dir: impl AsRef<Path>,
    ) -> Self {
        let step = step.into();
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis();
        let path = log_dir.as_ref().join(format!("{step}.{ts}.log"));
        Self::new_inner(spec, Some(ProgressLog { path, step }))
    }

    fn new_inner(spec: &'static StageSpec, log: Option<ProgressLog>) -> Self {
        Self {
            spec,
            counts: spec.metrics.iter().map(|_| AtomicUsize::new(0)).collect(),
            error_msg: Mutex::new(None),
            log,
            renderer: (!is_quiet()).then(Renderer::new),
        }
    }

    #[inline(always)]
    pub(crate) fn add(&self, idx: usize, count: usize) {
        self.counts[idx].fetch_add(count, Ordering::Relaxed);
    }

    #[inline(always)]
    pub(crate) fn inc(&self, idx: usize) {
        self.add(idx, 1);
    }

    /// `Annotating  1250000 reads | kept 1200000 | dropped 50000`
    fn line(&self, done: bool) -> String {
        let mut line = format!("{:<NAME_WIDTH$}", self.spec.name);
        if done {
            line.push_str("done: ");
        }
        line.push_str(&format!(
            "{} {}",
            self.counts[0].load(Ordering::Relaxed),
            self.spec.unit
        ));
        for (count, metric) in self.counts.iter().zip(self.spec.metrics.iter()).skip(1) {
            line.push_str(&format!(" | {metric} {}", count.load(Ordering::Relaxed)));
        }
        line
    }

    /// Cheap enough to call per read: it costs a clock read and an atomic compare
    /// until the redraw interval is actually due.
    #[inline(always)]
    pub(crate) fn refresh(&self) {
        if let Some(renderer) = &self.renderer
            && renderer.due()
        {
            renderer.draw(&self.line(false), false);
        }
    }

    pub(crate) fn store_error(&self, msg: impl Into<String>) {
        let msg = msg.into();
        self.print_error(msg.clone());
        if let Ok(mut err) = self.error_msg.lock() {
            *err = Some(msg);
        }
    }

    /// Errors are reported even under `--quiet`.
    pub(crate) fn print_error(&self, msg: impl Into<String>) {
        if let Some(renderer) = &self.renderer {
            renderer.erase();
        }
        eprintln!("{}", msg.into());
    }

    pub(crate) fn take_error(&self) -> Option<String> {
        self.error_msg.lock().ok().and_then(|mut err| err.take())
    }

    pub(crate) fn clear(&self) {
        if let Some(renderer) = &self.renderer {
            renderer.erase();
        }
    }

    pub(crate) fn finish(&self) {
        if let Some(log) = &self.log
            && let Err(e) = log.write(&self.counts, self.spec)
        {
            self.print_error(e);
        }

        if let Some(renderer) = &self.renderer {
            renderer.draw(&self.line(true), true);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_label_matches_legacy_format() {
        assert_eq!(log_label("total"), "Total:");
        assert_eq!(log_label("kept split"), "Kept split:");
        assert_eq!(log_label("failed"), "Failed:");
    }

    #[test]
    fn test_line_layout() {
        let tracker = ProgressTracker::new(&TRIM_STAGE);
        tracker.add(0, 100);
        tracker.add(1, 97);
        tracker.add(2, 4);
        tracker.add(3, 3);
        assert_eq!(
            tracker.line(false),
            "Trimming    100 reads | kept 97 | kept split 4 | failed 3"
        );
        assert_eq!(
            tracker.line(true),
            "Trimming    done: 100 reads | kept 97 | kept split 4 | failed 3"
        );
    }

    #[test]
    fn test_quiet_disables_rendering() {
        set_quiet(true);
        let tracker = ProgressTracker::new(&QC_STAGE);
        assert!(tracker.renderer.is_none());
        set_quiet(false);
        let tracker = ProgressTracker::new(&QC_STAGE);
        assert!(tracker.renderer.is_some());
    }
}
