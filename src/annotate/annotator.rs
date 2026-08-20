use crate::annotate::barcodes::{BarcodeGroup, BarcodeType};
use crate::annotate::edit_model::get_edit_cut_off;
use crate::annotate::searcher::{SarraceniaMatch, Demuxer};
use crate::config::AnnotateConfig;
use crate::io::io::open_fastq;
use crate::progress::progress::{ANNOTATE_STAGE, ProgressTracker, info, is_quiet};
use crate::qc::qc::{QcFilter, print_qc_summary, reset_qc_counts};
use anyhow::anyhow;
use seq_io::fastq::{Error as FastqError, Record, RecordSet};
use seq_io::parallel::{ParallelRecordsets, read_parallel};
use std::fmt::Display;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread_local;
use rayon::prelude::*;
use std::fs::File;

fn wrap_error<T, E: Display>(result: Result<T, E>, context: &str) -> anyhow::Result<T> {
    result.map_err(|err| anyhow!("{context}: {err}"))
}

#[inline(always)]
fn write_annotation_batch(
    writer: &Arc<Mutex<csv::Writer<std::fs::File>>>,
    record_set_results: &[SarraceniaMatch],
) -> anyhow::Result<()> {
    let mut writer = wrap_error(writer.lock(), "Annotation writer lock failed")?;
    for annotation in record_set_results {
        wrap_error(
            writer.serialize(annotation),
            "Failed to write annotation output",
        )?;
    }
    Ok(())
}

#[inline(always)]
fn consume_record_sets(
    record_sets: &mut ParallelRecordsets<RecordSet, FastqError, (Vec<SarraceniaMatch>, usize)>,
    writer: &Arc<Mutex<csv::Writer<std::fs::File>>>,
    progress: &ProgressTracker,
) -> anyhow::Result<()> {
    while let Some(result) = record_sets.next() {
        let (record_set, (record_set_results, found_count)) =
            wrap_error(result, "Input FASTQ parsing failed")?;
        let n_records = record_set.len();
        write_annotation_batch(writer, &record_set_results)?;
        progress.add(0, n_records);
        progress.add(1, found_count);
        progress.add(2, n_records - found_count);
        progress.refresh();
    }
    Ok(())
}

fn process_bam_batch(
    records: &[noodles::bam::Record],
    query_groups: &[BarcodeGroup],
    writer: &Arc<Mutex<csv::Writer<std::fs::File>>>,
    progress: &ProgressTracker,
    alpha: f32,
    verbose: bool,
    min_score: f64,
    min_score_diff: f64,
    qc: &QcFilter,
) -> anyhow::Result<()> {
    let results: Vec<(Vec<SarraceniaMatch>, bool)> = records
        .par_iter()
        .map(|record| {
            let read_id = match record.name() {
                Some(name) => String::from_utf8_lossy(name).into_owned(),
                None => return (Vec::new(), false),
            };

            let seq = record.sequence();

            // Pre-alignment gate: skip the demuxer entirely for reads that fail QC.
            if qc.is_active() {
                let quality = qc.quality_of_phred(record.quality_scores().iter());
                if !qc.passes(seq.len(), quality) {
                    return (Vec::new(), false);
                }
            }

            let mut seq_bytes = Vec::with_capacity(seq.len());
            for base in seq.iter() {
                seq_bytes.push(base.to_ascii_uppercase());
            }

            let mut demuxer = Demuxer::new(alpha, verbose, min_score, min_score_diff);
            for query_group in query_groups {
                demuxer.add_query_group(query_group.clone());
            }

            let matches = demuxer.demux(&read_id, &seq_bytes);
            let found = !matches.is_empty();
            (matches, found)
        })
        .collect();

    let mut all_matches = Vec::new();
    let mut found_count = 0;
    for (matches, found) in results {
        if found {
            found_count += 1;
            all_matches.extend(matches);
        }
    }

    write_annotation_batch(writer, &all_matches)?;
    let n_records = records.len();
    progress.add(0, n_records);
    progress.add(1, found_count);
    progress.add(2, n_records - found_count);
    progress.refresh();

    Ok(())
}

fn annotate_bam(
    read_file: &str,
    out_file: &str,
    query_groups: Vec<BarcodeGroup>,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    let alpha = config.alpha;
    let verbose = config.verbose;
    let min_score = config.min_score;
    let min_score_diff = config.min_score_diff;

    let mut reader = File::open(read_file)
        .map(noodles::bam::io::Reader::new)
        .map_err(|e| anyhow!("Failed to open BAM file '{read_file}': {e}"))?;

    let _header = reader.read_header()?;

    let writer = Arc::new(Mutex::new(
        csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_file)
            .map_err(|e| anyhow!("Failed to create annotation output file '{out_file}': {e}"))?,
    ));

    if !is_quiet() {
        for (i, query_group) in query_groups.iter().enumerate() {
            info(format!("{}: {}", query_group.barcode_type.as_str(), i));
            query_group.display(5);
        }
    }

    let progress = if config.verbose {
        let log_dir = Path::new(out_file)
            .parent()
            .unwrap_or_else(|| Path::new("."));
        ProgressTracker::new_with_logging(&ANNOTATE_STAGE, "annotate", log_dir)
    } else {
        ProgressTracker::new(&ANNOTATE_STAGE)
    };

    let mut record_buffer = Vec::with_capacity(1000);
    for result in reader.records() {
        let record = result.map_err(|e| anyhow!("Failed to read BAM record: {e}"))?;
        record_buffer.push(record);
        if record_buffer.len() >= 1000 {
            process_bam_batch(&record_buffer, &query_groups, &writer, &progress, alpha, verbose, min_score, min_score_diff, &config.qc)?;
            record_buffer.clear();
        }
    }
    if !record_buffer.is_empty() {
        process_bam_batch(&record_buffer, &query_groups, &writer, &progress, alpha, verbose, min_score, min_score_diff, &config.qc)?;
    }

    progress.finish();
    print_qc_summary(&config.qc);
    Ok(())
}

// used by custom experiments (direct annotate call)
pub fn annotate_with_files(
    read_file: &str,
    query_files: Vec<&str>,
    query_types: Vec<BarcodeType>,
    out_file: &str,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    // Get query groups
    let mut query_groups = Vec::new();
    for (query_file, query_type) in query_files.iter().zip(query_types.iter()) {
        let mut query_group = BarcodeGroup::new_from_fasta(query_file, query_type.clone());
        if let Some(max_flank_errors) = config.max_flank_errors {
            query_group.set_flank_threshold(max_flank_errors);
        } else {
            // Determine based on formula
            let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
            info(format!("Auto edit flank cut off: {edit_cut_off}"));
            query_group.set_flank_threshold(edit_cut_off);
        }
        query_groups.push(query_group);
    }
    annotate(read_file, out_file, query_groups, config)
}

// we could maybe just discard annotate_with_groups and only have kit or fasta
pub fn annotate_with_kit(
    read_file: &str,
    out_file: &str,
    kit: &str,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    let query_groups: Vec<BarcodeGroup> = BarcodeGroup::new_from_kit(kit, config.use_extended);
    annotate_with_groups(read_file, out_file, query_groups, config)
}

// used by kit
pub fn annotate_with_groups(
    read_file: &str,
    out_file: &str,
    query_groups: Vec<BarcodeGroup>,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    let query_groups: Vec<BarcodeGroup> = query_groups
        .into_iter()
        .map(|mut query_group| {
            if let Some(max_flank_errors) = config.max_flank_errors {
                query_group.set_flank_threshold(max_flank_errors);
            } else {
                // Determine based on formula
                let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
                info(format!("Auto edit flank cut off: {edit_cut_off}"));
                query_group.set_flank_threshold(edit_cut_off);
            }
            query_group
        })
        .collect();
    annotate(read_file, out_file, query_groups, config)
}

pub fn annotate(
    read_file: &str,
    out_file: &str,
    query_groups: Vec<BarcodeGroup>,
    config: &AnnotateConfig,
) -> anyhow::Result<()> {
    config.qc.validate()?;
    reset_qc_counts();

    if read_file.to_ascii_lowercase().ends_with(".bam") {
        return annotate_bam(read_file, out_file, query_groups, config);
    }

    let alpha = config.alpha;
    let n_threads = config.n_threads;
    let verbose = config.verbose;
    let min_score = config.min_score;
    let min_score_diff = config.min_score_diff;
    let qc = config.qc;

    let reader = open_fastq(read_file);
    let writer = Arc::new(Mutex::new(
        csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(out_file)
            .map_err(|e| anyhow!("Failed to create annotation output file '{out_file}': {e}"))?,
    ));

    // Dispaly to user
    if !is_quiet() {
        for (i, query_group) in query_groups.iter().enumerate() {
            info(format!("{}: {}", query_group.barcode_type.as_str(), i));
            query_group.display(5);
        }
    }

    let progress = if config.verbose {
        let log_dir = Path::new(out_file)
            .parent()
            .unwrap_or_else(|| Path::new("."));
        ProgressTracker::new_with_logging(&ANNOTATE_STAGE, "annotate", log_dir)
    } else {
        ProgressTracker::new(&ANNOTATE_STAGE)
    };

    read_parallel(
        reader,
        n_threads,
        1000,
        |record_set| {
            // Create thread local demuxer if not init for current thread yet
            thread_local! {
                static DEMUXER: std::cell::RefCell<Option<Demuxer>> = const { std::cell::RefCell::new(None) };
            }
            DEMUXER.with(|cell| {
                if cell.borrow().is_none() {
                    let mut demux = Demuxer::new(alpha, verbose, min_score, min_score_diff);
                    for query_group in query_groups.iter() {
                        demux.add_query_group(query_group.clone());
                    }

                    *cell.borrow_mut() = Some(demux);
                }
            });

            // Go over the
            let mut record_set_annotations = Vec::new();
            let mut found = 0;
            for record in record_set.into_iter() {
                // Pre-alignment gate: skip the demuxer entirely for reads that fail QC.
                if qc.is_active() {
                    let quality = qc.quality_of_ascii(record.qual());
                    if !qc.passes(record.seq().len(), quality) {
                        continue;
                    }
                }

                // Use the demuxer through thread-local storage
                let matches: Vec<SarraceniaMatch> = DEMUXER.with(|cell| {
                    if let Some(ref mut demuxer) = *cell.borrow_mut() {
                        match record.id() {
                            Ok(read_id) => demuxer.demux(read_id, record.seq()),
                            Err(_) => Vec::new(),
                        }
                    } else {
                        Vec::new()
                    }
                });

                if !matches.is_empty() {
                    found += 1;
                    record_set_annotations.extend(matches);
                }
            }
            (record_set_annotations, found)
        },
        |record_sets| {
            if let Err(e) = consume_record_sets(record_sets, &writer, &progress) {
                progress.store_error(e.to_string());
            }
        },
    );

    if let Some(msg) = progress.take_error() {
        progress.clear();
        return Err(anyhow!(msg));
    }

    progress.finish();
    print_qc_summary(&qc);

    Ok(())
}
