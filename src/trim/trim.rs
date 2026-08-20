use crate::annotate::barcodes::BarcodeType;
use crate::annotate::searcher::SarraceniaMatch;
use crate::config::TrimConfig;
use crate::filter::pattern::{Cut, CutDirection};
use crate::io::io::open_fastq;
use crate::progress::progress::{ProgressTracker, TRIM_STAGE};
use crate::qc::qc::{print_qc_summary_labeled, reset_qc_counts};
use anyhow::anyhow;
use csv;
use flate2::Compression;
use flate2::write::GzEncoder;
use sassy::Strand;
use seq_io::fastq::Record;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use clap::ValueEnum;

use crate::trim::mod_recalculator::recalculate_base_mods;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::Header;
use noodles::sam::alignment::io::Write as SamWrite;


const TOTAL_IDX: usize = 0;
const TRIMMED_IDX: usize = 1;
const TRIMMED_SPLIT_IDX: usize = 2;
const FAILED_IDX: usize = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum LabelSide {
    Left,
    Right,
}

pub struct LabelConfig {
    pub include_label: bool,
    pub include_orientation: bool,
    pub include_flank: bool,
    pub sort_labels: bool,
    pub only_side: Option<LabelSide>,
}

impl LabelConfig {
    pub fn new(
        include_label: bool,
        include_orientation: bool,
        include_flank: bool,
        sort_labels: bool,
        only_side: Option<LabelSide>,
    ) -> Self {
        Self {
            include_label,
            include_orientation,
            include_flank,
            sort_labels,
            only_side,
        }
    }

    pub fn create_label(&self, annotations: &[SarraceniaMatch]) -> String {
        if !self.include_label {
            return "none".to_string();
        }

        let mut label_parts: Vec<String> = annotations
            .iter()
            .filter_map(|m| {
                let label = m.label.clone();

                // Skip if it's a flank and we don't want flanks
                // this also prevents having a flank come before label when
                // selecting only left or only right
                if !self.include_flank && label.contains("flank") {
                    return None;
                }

                let mut result = label;

                if self.include_orientation {
                    let ori = match m.strand {
                        Strand::Fwd => "fw",
                        Strand::Rc => "rc",
                    };
                    result = format!("{result}_{ori}");
                }
                Some(result)
            })
            .collect();

        if self.sort_labels && self.only_side.is_some() {
            panic!("Cannot enable only keeping left label and sorting as this makes it ambiguous");
        }

        if label_parts.is_empty() {
            "none".to_string()
        } else if self.sort_labels {
            label_parts.sort();
            label_parts.join("__")
        } else if self.only_side.is_some() {
            let side = self.only_side.unwrap();
            if side == LabelSide::Left {
                label_parts.first().unwrap().clone()
            } else {
                label_parts.last().unwrap().clone()
            }
        } else {
            label_parts.join("__")
        }
    }
}

impl TrimConfig {
    pub fn get_label_config(&self) -> LabelConfig {
        LabelConfig::new(
            self.add_labels,
            self.add_orientation,
            self.add_flank,
            self.sort_labels,
            self.only_side,
        )
    }
}

#[derive(Debug)]
struct CompleteSlice {
    start: usize,
    end: usize,
    annotations: Vec<SarraceniaMatch>,
}

fn preprocess_cuts(annotations: &[SarraceniaMatch], seq_len: usize) -> Vec<CompleteSlice> {
    let mut slices: Vec<CompleteSlice> = Vec::new();

    // Group cuts by their IDs
    let mut cut_groups: HashMap<usize, Vec<(usize, usize, &Cut, &SarraceniaMatch)>> = HashMap::new();
    for anno in annotations {
        if let Some(cuts) = &anno.cuts {
            for (cut, _) in cuts {
                cut_groups.entry(cut.group_id).or_default().push((
                    anno.read_start_flank,
                    anno.read_end_flank,
                    cut,
                    anno,
                ));
            }
        }
    }

    // Sort groups by their leftmost position
    let mut sorted_groups: Vec<_> = cut_groups.into_iter().collect();
    sorted_groups.sort_by_key(|(_, group)| {
        group
            .first()
            .map(|(start, _, _, _)| *start)
            .unwrap_or(usize::MAX)
    });

    // Process each group
    for (i, (_, group)) in sorted_groups.iter().enumerate() {
        if group.len() == 2 {
            // We have two annotations so get start and end based on their cuts
            let group1 = &group[0];
            let group2 = &group[1];

            // Get start position based on first group's cut direction
            let start = match &group1.2.direction {
                CutDirection::Before => group1.0,
                CutDirection::After => group1.1,
            };

            // Get end position based on second group's cut direction
            let end = match &group2.2.direction {
                CutDirection::Before => group2.0,
                CutDirection::After => group2.1,
            };

            let annotations = vec![group1.3.clone(), group2.3.clone()];

            slices.push(CompleteSlice {
                start,
                end,
                annotations,
            });
        } else if group.len() == 1 {
            let &(start, end, cut, anno) = &group[0];

            match cut.direction {
                CutDirection::Before => {
                    // Look left for start position and annotation
                    let (slice_start, left_anno) = if i > 0 {
                        let prev_group = &sorted_groups[i - 1].1;
                        let max_end_idx = prev_group
                            .iter()
                            .enumerate()
                            .max_by_key(|(_, (_, end, _, _))| end)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (
                            prev_group[max_end_idx].1,
                            Some(prev_group[max_end_idx].3.clone()),
                        )
                    } else {
                        (0, None)
                    };

                    let mut annotations = Vec::new();
                    if let Some(left) = left_anno {
                        annotations.push(left);
                    }
                    annotations.push(anno.clone());

                    slices.push(CompleteSlice {
                        start: slice_start,
                        end: start,
                        annotations,
                    });
                }
                CutDirection::After => {
                    // Look right for end position and annotation
                    let (slice_end, right_anno) = if i < sorted_groups.len() - 1 {
                        let next_group = &sorted_groups[i + 1].1;
                        let min_start_idx = next_group
                            .iter()
                            .enumerate()
                            .min_by_key(|(_, (start, _, _, _))| start)
                            .map(|(idx, _)| idx)
                            .unwrap_or(0);
                        (
                            next_group[min_start_idx].0,
                            Some(next_group[min_start_idx].3.clone()),
                        )
                    } else {
                        (seq_len, None)
                    };

                    let mut annotations = Vec::new();
                    annotations.push(anno.clone());
                    if let Some(right) = right_anno {
                        annotations.push(right);
                    }

                    slices.push(CompleteSlice {
                        start: end,
                        end: slice_end,
                        annotations,
                    });
                }
            }
        }
    }
    slices
}

pub fn process_read_and_anno(
    seq: &[u8],
    qual: &[u8],
    annotations: &[SarraceniaMatch],
    label_config: &LabelConfig,
    skip_trim: bool,
    flip: bool,
) -> Vec<(Vec<u8>, Vec<u8>, String, String)> {
    let mut results = Vec::new();
    let seq_len = seq.len();

    // Preprocess cuts to get complete slices
    let slices = preprocess_cuts(annotations, seq_len);

    // Group slices by cut group ID
    for (slice_count, slice) in slices.iter().enumerate() {
        if slice.start >= slice.end {
            continue;
        }

        // For now if trimming is disabled, we just
        // return the full sequence and quality
        let mut trimmed_seq = if skip_trim {
            seq.to_vec()
        } else {
            seq[slice.start..slice.end].to_vec()
        };
        let mut trimmed_qual = if skip_trim {
            qual.to_vec()
        } else {
            qual[slice.start..slice.end].to_vec()
        };

        if flip && should_flip(&slice.annotations) {
            trimmed_seq = reverse_complement(&trimmed_seq);
            trimmed_qual.reverse();
        }

        let label_matches: Vec<SarraceniaMatch> = slice.annotations.clone();

        let group_label = label_config.create_label(&label_matches);
        let read_suffix = if slice_count == 0 {
            "".to_string()
        } else {
            format!("_{slice_count}")
        };
        results.push((trimmed_seq, trimmed_qual, group_label, read_suffix));
    }

    results
}

/// Extracts the clean read ID from a FASTQ record ID by taking the first part before any whitespace
#[allow(unused)]
fn clean_read_id(id: &str) -> &str {
    id.split_whitespace().next().unwrap_or(id)
}

fn format_output_file_error(output_file: &str, err: &std::io::Error) -> String {
    let mut msg = format!("Failed to create output file '{output_file}': {err}");
    if err.raw_os_error() == Some(24) {
        msg.push_str("\nTry setting ulimit higher: \"ulimit -n 65000\"");
    }
    msg
}

/// Record a read that produced no output, either because trimming yielded nothing or
/// because everything it yielded was dropped by the post-trim QC gate.
fn mark_failed(
    progress: &ProgressTracker,
    failed_writer: &mut Option<BufWriter<File>>,
    read_id: &str,
) {
    progress.inc(FAILED_IDX);
    if let Some(writer) = failed_writer.as_mut() {
        writeln!(writer, "{read_id}").expect("Failed to write to failed trimmed writer");
    }
}

fn should_flip(annotations: &[SarraceniaMatch]) -> bool {
    // If we matched an Ftag in rc we flip
    annotations
        .iter()
        .any(|anno| anno.match_type == BarcodeType::Ftag && anno.strand == Strand::Rc)
}

pub fn trim_matches(
    filtered_match_file: &str,
    read_fastq_file: &str,
    output_folder: &str,
    config: &TrimConfig,
) -> anyhow::Result<()> {
    config.qc.validate_flags("trimmed-")?;
    if config.qc.is_active() {
        reset_qc_counts();
    }

    if read_fastq_file.to_ascii_lowercase().ends_with(".bam") {
        return trim_bam_matches(filtered_match_file, read_fastq_file, output_folder, config);
    }

    // Create output folder if it doesn't exist
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    // Label formatting config
    let label_config = config.get_label_config();

    if config.sort_labels && config.only_side.is_some() {
        return Err(anyhow!(
            "Cannot enable only keeping left/right label and sorting; this is ambiguous"
        ));
    }

    // Read all annotations and group by read ID
    let mut annotations_by_read: HashMap<String, Vec<SarraceniaMatch>> = HashMap::new();

    // Create progress bars
    let progress = if config.verbose {
        ProgressTracker::new_with_logging(&TRIM_STAGE, "trim", output_folder)
    } else {
        ProgressTracker::new(&TRIM_STAGE)
    };

    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    for result in matches_reader.deserialize() {
        let anno: SarraceniaMatch = result.expect("Failed to parse annotation line");
        annotations_by_read
            .entry(anno.read_id.clone())
            .or_default()
            .push(anno);
    }

    // Create writers regular or gzip write
    let mut writers: HashMap<String, Box<dyn Write>> = HashMap::new();

    // If there is a failed trimmed writer, create it
    let mut failed_trimmed_writer =
        config
            .failed_trimmed_writer
            .as_ref()
            .map(|failed_trimmed_writer_path| {
                BufWriter::new(File::create(failed_trimmed_writer_path).unwrap())
            });

    // Process reads
    let mut reader = open_fastq(read_fastq_file);

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");
        let (read_id, desc) = record.id_desc().unwrap();
        let read_id = read_id.to_string();
        let desc: &str = desc.unwrap_or_default();
        progress.inc(TOTAL_IDX);

        // Check if this read has annotations
        if let Some(annotations) = annotations_by_read.get(&read_id) {
            // mapped_reads += 1;

            let mut results: Vec<(Vec<u8>, Vec<u8>, String, String)> = process_read_and_anno(
                record.seq(),
                record.qual(),
                annotations,
                &label_config,
                config.skip_trim,
                config.flip,
            );

            // Post-trim gate: judge what we are about to write, not the raw read.
            if config.qc.is_active() {
                results.retain(|(seq, qual, _, _)| config.qc.passes_ascii(seq.len(), qual));
            }

            if !results.is_empty() {
                progress.inc(TRIMMED_IDX);
            } else {
                mark_failed(&progress, &mut failed_trimmed_writer, &read_id);
            }

            if results.len() > 1 {
                progress.inc(TRIMMED_SPLIT_IDX);
            }

            for (trimmed_seq, trimmed_qual, group, read_suffix) in results {
                // Get or create writer for this group
                if !writers.contains_key(&group) {
                    let output_file = if config.gzip {
                        format!("{output_folder}/{group}.trimmed.fastq.gz")
                    } else {
                        format!("{output_folder}/{group}.trimmed.fastq")
                    };
                    let file = File::create(&output_file).map_err(|err| {
                        let msg = format_output_file_error(&output_file, &err);
                        progress.print_error(msg.clone());
                        progress.clear();
                        anyhow!(msg)
                    })?;

                    // If gzipping is enabled we use zipped writer
                    let writer: Box<dyn Write> = if config.gzip {
                        Box::new(GzEncoder::new(file, Compression::default()))
                    } else {
                        Box::new(BufWriter::new(file))
                    };
                    writers.insert(group.clone(), writer);
                }
                let writer = writers
                    .get_mut(&group)
                    .expect("writer should exist after insertion");

                // Write FASTQ format
                if config.write_full_header {
                    writeln!(writer, "@{read_id}{read_suffix} {desc}")
                        .expect("Failed to write header");
                } else {
                    writeln!(writer, "@{read_id}{read_suffix}").expect("Failed to write header");
                }
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_seq))
                    .expect("Failed to write sequence");
                writeln!(writer, "+").expect("Failed to write separator");
                writeln!(writer, "{}", String::from_utf8_lossy(&trimmed_qual))
                    .expect("Failed to write quality");
            }
        }

        progress.refresh();
    }

    // Flush all writers
    for (_, writer) in writers.iter_mut() {
        writer.flush().expect("Failed to flush output");
    }

    progress.finish();
    if config.qc.is_active() {
        print_qc_summary_labeled(&config.qc, "Post-trim QC");
    }
    Ok(())
}

#[inline(always)]
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| RC[c as usize]).collect()
}

const RC: [u8; 256] = {
    let mut rc = [0; 256];
    let mut i = 0;
    while i < 256 {
        rc[i] = i as u8;
        i += 1;
    }
    // Standard bases
    rc[b'A' as usize] = b'T';
    rc[b'C' as usize] = b'G';
    rc[b'T' as usize] = b'A';
    rc[b'G' as usize] = b'C';
    rc[b'a' as usize] = b't';
    rc[b'c' as usize] = b'g';
    rc[b't' as usize] = b'a';
    rc[b'g' as usize] = b'c';
    // IUPAC ambiguity codes
    rc[b'R' as usize] = b'Y'; // A|G -> T|C
    rc[b'Y' as usize] = b'R'; // C|T -> G|A
    rc[b'S' as usize] = b'S'; // G|C -> C|G
    rc[b'W' as usize] = b'W'; // A|T -> T|A
    rc[b'K' as usize] = b'M'; // G|T -> C|A
    rc[b'M' as usize] = b'K'; // A|C -> T|G
    rc[b'B' as usize] = b'V'; // C|G|T -> G|C|A
    rc[b'D' as usize] = b'H'; // A|G|T -> T|C|A
    rc[b'H' as usize] = b'D'; // A|C|T -> T|G|A
    rc[b'V' as usize] = b'B'; // A|C|G -> T|G|C
    rc[b'N' as usize] = b'N'; // A|C|G|T -> T|G|C|A
    rc[b'X' as usize] = b'X';
    // Lowercase versions
    rc[b'r' as usize] = b'y';
    rc[b'y' as usize] = b'r';
    rc[b's' as usize] = b's';
    rc[b'w' as usize] = b'w';
    rc[b'k' as usize] = b'm';
    rc[b'm' as usize] = b'k';
    rc[b'b' as usize] = b'v';
    rc[b'd' as usize] = b'h';
    rc[b'h' as usize] = b'd';
    rc[b'v' as usize] = b'b';
    rc[b'n' as usize] = b'n';
    rc[b'x' as usize] = b'x';
    rc
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotate::barcodes::BarcodeType;
    use crate::filter::pattern::{Cut, CutDirection};

    #[test]
    fn test_single_cut() {
        let seq = b"CCCCCCCCAAAACCCCCCCCCCCC";
        let qual = b"________IIII____________";

        let annotations = vec![
            SarraceniaMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "Fbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::After), 8)]),
            ),
            SarraceniaMatch::new(
                12, // read_start_bar
                16, // read_end_bar
                12, // read_start_flank
                16, // read_end_flank
                0,  // bar_start
                4,  // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "Rbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::Before), 12)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AAAA");
        assert_eq!(trimmed_qual, b"IIII");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");
    }

    #[test]
    fn test_two_cut_groups_produce_two_slices() {
        // seq indices: 0..8 C, 8..20 A, 20..26 C, 26..28 G, 28..30 C
        let seq = b"CCCCCCCCAAAAAAAAAAAACCCCCCGGCC";
        let qual = b"________IIIIIIIIIIII______II__";

        let read_len = seq.len();

        let annotations = vec![
            // Group 1: start at After(end_flank=8), end at Before(start_flank=20) -> slice 8..20
            SarraceniaMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "F1".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(1, CutDirection::After), 8)]),
            ),
            SarraceniaMatch::new(
                20, // read_start_bar
                24, // read_end_bar
                20, // read_start_flank
                24, // read_end_flank
                0,
                4,
                BarcodeType::Rtag,
                0,
                0,
                "R1".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(1, CutDirection::Before), 20)]),
            ),
            // Group 2: start at After(end_flank=26), end at Before(start_flank=28) -> slice 26..28
            SarraceniaMatch::new(
                24,
                26,
                24,
                26,
                0,
                2,
                BarcodeType::Ftag,
                0,
                0,
                "F2".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(2, CutDirection::After), 26)]),
            ),
            SarraceniaMatch::new(
                28,
                30,
                28,
                30,
                0,
                2,
                BarcodeType::Rtag,
                0,
                0,
                "R2".to_string(),
                Strand::Fwd,
                read_len,
                "read1".to_string(),
                0,
                Some(vec![(Cut::new(2, CutDirection::Before), 28)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, false);

        assert_eq!(results.len(), 2);

        let (trimmed_seq1, trimmed_qual1, label1, _) = &results[0];
        assert_eq!(trimmed_seq1, b"AAAAAAAAAAAA");
        assert_eq!(trimmed_qual1, b"IIIIIIIIIIII");
        assert_eq!(label1, "F1_fw__R1_fw");

        let (trimmed_seq2, trimmed_qual2, label2, _) = &results[1];
        assert_eq!(trimmed_seq2, b"GG");
        assert_eq!(trimmed_qual2, b"II");
        assert_eq!(label2, "F2_fw__R2_fw");
    }

    #[test]
    fn test_trim_skipping() {
        let seq = b"CCCCCCCCAAAACCCCCCCCCCCC";
        let qual = b"________IIII____________";

        let annotations = vec![
            SarraceniaMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "Fbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::After), 8)]),
            ),
            SarraceniaMatch::new(
                12, // read_start_bar
                16, // read_end_bar
                12, // read_start_flank
                16, // read_end_flank
                0,  // bar_start
                4,  // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "Rbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::Before), 12)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, true, false);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"CCCCCCCCAAAACCCCCCCCCCCC");
        assert_eq!(trimmed_qual, b"________IIII____________");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");
    }

    #[test]
    fn test_flipping() {
        let seq = b"CCCCCCCCAGGCCCCCCCCCCCCC";
        let qual = b"________IIIA____________";

        let mut annotations = vec![
            SarraceniaMatch::new(
                4, // read_start_bar
                8, // read_end_bar
                4, // read_start_flank
                8, // read_end_flank
                0, // bar_start
                4, // bar_end
                BarcodeType::Ftag,
                0, // flank_cost
                0, // barcode_cost
                "Fbar".to_string(),
                Strand::Rc, // Note RC match we have to flip
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::After), 8)]),
            ),
            SarraceniaMatch::new(
                12, // read_start_bar
                16, // read_end_bar
                12, // read_start_flank
                16, // read_end_flank
                0,  // bar_start
                4,  // bar_end
                BarcodeType::Rtag,
                0, // flank_cost
                0, // barcode_cost
                "Rbar".to_string(),
                Strand::Fwd,
                seq.len(),
                "read1".to_string(),
                0, // rel_dist_to_end
                Some(vec![(Cut::new(0, CutDirection::Before), 12)]),
            ),
        ];

        let label_config = LabelConfig::new(true, true, true, true, None);
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true);

        assert_eq!(results.len(), 1);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"GCCT");
        assert_eq!(trimmed_qual, b"AIII");
        assert_eq!(group_label, "Fbar_rc__Rbar_fw");

        annotations[0].strand = Strand::Fwd;
        let results = process_read_and_anno(seq, qual, &annotations, &label_config, false, true);
        let (trimmed_seq, trimmed_qual, group_label, _) = &results[0];
        println!("trimmed_seq: {}", String::from_utf8_lossy(trimmed_seq));
        assert_eq!(trimmed_seq, b"AGGC");
        assert_eq!(trimmed_qual, b"IIIA");
        assert_eq!(group_label, "Fbar_fw__Rbar_fw");

        // Chaning the strand in Fbar match should give original seq and qual again
    }

    use crate::qc::qc::{QcFilter, lock_qc_counts_for_test, qc_counts};
    use tempfile::TempDir;

    /// Two matches that cut `CCCCCCCCAAAACCCCCCCCCCCC` down to the 4bp `AAAA` insert.
    fn insert_annotations(read_id: &str, read_len: usize) -> Vec<SarraceniaMatch> {
        vec![
            SarraceniaMatch::new(
                4,
                8,
                4,
                8,
                0,
                4,
                BarcodeType::Ftag,
                0,
                0,
                "Fbar".to_string(),
                Strand::Fwd,
                read_len,
                read_id.to_string(),
                0,
                Some(vec![(Cut::new(0, CutDirection::After), 8)]),
            ),
            SarraceniaMatch::new(
                12,
                16,
                12,
                16,
                0,
                4,
                BarcodeType::Rtag,
                0,
                0,
                "Rbar".to_string(),
                Strand::Fwd,
                read_len,
                read_id.to_string(),
                0,
                Some(vec![(Cut::new(0, CutDirection::Before), 12)]),
            ),
        ]
    }

    fn trim_config_with_qc(qc: QcFilter, failed_out: Option<String>) -> TrimConfig {
        TrimConfig {
            add_labels: true,
            add_orientation: false,
            add_flank: true,
            sort_labels: false,
            only_side: None,
            failed_trimmed_writer: failed_out,
            write_full_header: true,
            skip_trim: false,
            flip: false,
            verbose: false,
            gzip: false,
            qc,
        }
    }

    /// Writes a two-read FASTQ (one high-quality insert, one low-quality insert) plus
    /// the matching annotation TSV, and returns `(dir, matches_path, reads_path)`.
    fn post_trim_fixture() -> (TempDir, String, String) {
        let dir = TempDir::new().unwrap();

        let reads_path = dir.path().join("reads.fastq");
        std::fs::write(
            &reads_path,
            b"@hq_read\nCCCCCCCCAAAACCCCCCCCCCCC\n+\n________IIII____________\n\
              @lq_read\nCCCCCCCCAAAACCCCCCCCCCCC\n+\n________!!!!____________\n"
                .as_slice(),
        )
        .unwrap();

        let matches_path = dir.path().join("filtered.tsv");
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(&matches_path)
            .unwrap();
        for read_id in ["hq_read", "lq_read"] {
            for anno in insert_annotations(read_id, 24) {
                writer.serialize(anno).unwrap();
            }
        }
        writer.flush().unwrap();

        let matches = matches_path.to_str().unwrap().to_string();
        let reads = reads_path.to_str().unwrap().to_string();
        (dir, matches, reads)
    }

    #[test]
    fn test_post_trim_quality_gate_drops_low_quality_insert() {
        let _guard = lock_qc_counts_for_test();
        let (dir, matches, reads) = post_trim_fixture();
        let out = dir.path().join("out");
        let failed = dir.path().join("failed.txt");

        // Both reads are 24bp raw and would sail through a pre-alignment gate; only the
        // 4bp trimmed inserts differ in quality.
        let config = trim_config_with_qc(
            QcFilter::new(None, None, Some(10.0)),
            Some(failed.to_str().unwrap().to_string()),
        );
        trim_matches(&matches, &reads, out.to_str().unwrap(), &config).unwrap();

        let trimmed = std::fs::read_to_string(out.join("Fbar__Rbar.trimmed.fastq")).unwrap();
        assert!(trimmed.contains("@hq_read"), "got: {trimmed}");
        assert!(!trimmed.contains("@lq_read"), "got: {trimmed}");

        // A read whose only fragment is rejected counts as failed, not trimmed.
        let failed_ids = std::fs::read_to_string(&failed).unwrap();
        assert_eq!(failed_ids.trim(), "lq_read");
        assert_eq!(qc_counts().dropped_quality, 1);
    }

    #[test]
    fn test_post_trim_length_gate_uses_trimmed_length() {
        let _guard = lock_qc_counts_for_test();
        let (dir, matches, reads) = post_trim_fixture();
        let out = dir.path().join("out");

        // 5 > the 4bp insert but well under the 24bp raw read: a pre-alignment gate
        // would have kept both reads.
        let config = trim_config_with_qc(QcFilter::new(Some(5), None, None), None);
        trim_matches(&matches, &reads, out.to_str().unwrap(), &config).unwrap();

        // Writers are opened lazily, so rejecting everything leaves no output file at all.
        assert!(!out.join("Fbar__Rbar.trimmed.fastq").exists());
        assert_eq!(qc_counts().dropped_length, 2);
    }

    #[test]
    fn test_post_trim_gate_off_keeps_everything() {
        let _guard = lock_qc_counts_for_test();
        let (dir, matches, reads) = post_trim_fixture();
        let out = dir.path().join("out");

        let config = trim_config_with_qc(QcFilter::default(), None);
        trim_matches(&matches, &reads, out.to_str().unwrap(), &config).unwrap();

        let trimmed = std::fs::read_to_string(out.join("Fbar__Rbar.trimmed.fastq")).unwrap();
        assert!(trimmed.contains("@hq_read"));
        assert!(trimmed.contains("@lq_read"));
        assert_eq!(qc_counts().dropped_length + qc_counts().dropped_quality, 0);
    }

    #[test]
    fn test_post_trim_rejects_impossible_bounds() {
        let (dir, matches, reads) = post_trim_fixture();
        let out = dir.path().join("out");
        let config = trim_config_with_qc(QcFilter::new(Some(500), Some(100), None), None);
        assert!(trim_matches(&matches, &reads, out.to_str().unwrap(), &config).is_err());
    }
}

pub fn trim_bam_matches(
    filtered_match_file: &str,
    read_bam_file: &str,
    output_folder: &str,
    config: &TrimConfig,
) -> anyhow::Result<()> {
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
    }

    let label_config = config.get_label_config();

    // Read all annotations and group by read ID
    let mut annotations_by_read: HashMap<String, Vec<SarraceniaMatch>> = HashMap::new();

    let progress = if config.verbose {
        ProgressTracker::new_with_logging(&TRIM_STAGE, "trim", output_folder)
    } else {
        ProgressTracker::new(&TRIM_STAGE)
    };

    let mut matches_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(filtered_match_file)
        .expect("Failed to open matches file");

    for result in matches_reader.deserialize() {
        let anno: SarraceniaMatch = result.expect("Failed to parse annotation line");
        annotations_by_read
            .entry(anno.read_id.clone())
            .or_default()
            .push(anno);
    }

    let mut file_reader = File::open(read_bam_file)
        .map(noodles::bam::io::Reader::new)
        .map_err(|e| anyhow!("Failed to open input BAM: {e}"))?;

    let header = file_reader.read_header()?;

    // Output writers map (using precise BGZF wrapped BAM writer type)
    let mut writers: HashMap<String, noodles::bam::io::Writer<noodles::bgzf::io::Writer<File>>> = HashMap::new();

    let mut failed_trimmed_writer =
        config
            .failed_trimmed_writer
            .as_ref()
            .map(|failed_trimmed_writer_path| {
                BufWriter::new(File::create(failed_trimmed_writer_path).unwrap())
            });

    let mut warned_about_alignment = false;

    for result in file_reader.records() {
        let record = result.map_err(|e| anyhow!("Error reading record: {e}"))?;
        let read_id = match record.name() {
            Some(name) => String::from_utf8_lossy(name).into_owned(),
            None => continue,
        };
        progress.inc(TOTAL_IDX);

        if let Some(annotations) = annotations_by_read.get(&read_id) {
            let seq = record.sequence();
            let seq_bytes: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();

            let qual = record.quality_scores();
            let qual_bytes: Vec<u8> = qual.iter().collect();

            // Get original MM/ML tags
            let mm_tag = Tag::try_from([b'M', b'M']).unwrap();
            let ml_tag = Tag::try_from([b'M', b'L']).unwrap();

            let orig_mm = record.data().get(&mm_tag).and_then(|v_res| {
                if let Ok(noodles::sam::alignment::record::data::field::Value::String(ref s)) = v_res {
                    std::str::from_utf8(s.as_ref()).ok().map(|str_slice| str_slice.to_string())
                } else if let Ok(noodles::sam::alignment::record::data::field::Value::Character(c)) = v_res {
                    Some((c as char).to_string())
                } else {
                    None
                }
            });

            let orig_ml = record.data().get(&ml_tag).and_then(|v_res| {
                if let Ok(noodles::sam::alignment::record::data::field::Value::Array(noodles::sam::alignment::record::data::field::value::Array::UInt8(ref arr))) = v_res {
                    let mut probs = Vec::new();
                    for val_res in arr.iter() {
                        if let Ok(val) = val_res {
                            probs.push(val);
                        }
                    }
                    Some(probs)
                } else if let Ok(noodles::sam::alignment::record::data::field::Value::Array(noodles::sam::alignment::record::data::field::value::Array::Int8(ref arr))) = v_res {
                    let mut probs = Vec::new();
                    for val_res in arr.iter() {
                        if let Ok(val) = val_res {
                            probs.push(val as u8);
                        }
                    }
                    Some(probs)
                } else {
                    None
                }
            });

            // Find slices/cuts
            let slices = preprocess_cuts(annotations, seq_bytes.len());

            // Counted rather than derived from `slices`, because a slice can still be
            // discarded below (empty span, or rejected by the post-trim QC gate).
            let mut written = 0usize;

            for (slice_count, slice) in slices.iter().enumerate() {
                if slice.start >= slice.end {
                    continue;
                }

                let mut trimmed_seq = if config.skip_trim {
                    seq_bytes.clone()
                } else {
                    seq_bytes[slice.start..slice.end].to_vec()
                };

                let mut trimmed_qual = if config.skip_trim {
                    qual_bytes.clone()
                } else {
                    qual_bytes[slice.start..slice.end].to_vec()
                };

                let flip_read = config.flip && should_flip(&slice.annotations);
                if flip_read {
                    trimmed_seq = reverse_complement(&trimmed_seq);
                    trimmed_qual.reverse();
                }

                // Post-trim gate: judge what we are about to write, not the raw read.
                // Checked before the MM/ML recalculation so rejects cost nothing.
                if config.qc.is_active()
                    && !config.qc.passes_phred(trimmed_seq.len(), &trimmed_qual)
                {
                    continue;
                }

                // Recalculate MM/ML
                let mut new_mm = None;
                let mut new_ml = None;
                if let (Some(mm), Some(ml)) = (&orig_mm, &orig_ml) {
                    let (mm_recalc, ml_recalc) = recalculate_base_mods(
                        &seq_bytes,
                        &trimmed_seq,
                        if config.skip_trim { 0 } else { slice.start },
                        if config.skip_trim { seq_bytes.len() } else { slice.end },
                        flip_read,
                        mm,
                        ml,
                    );
                    if !mm_recalc.is_empty() {
                        new_mm = Some(mm_recalc);
                        new_ml = Some(ml_recalc);
                    }
                }

                // Create trimmed record buffer
                let mut record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;

                // Update sequence & quality
                let seq_buf = noodles::sam::alignment::record_buf::Sequence::from(trimmed_seq.clone());
                let qual_buf = noodles::sam::alignment::record_buf::QualityScores::from(trimmed_qual);
                *record_buf.sequence_mut() = seq_buf;
                *record_buf.quality_scores_mut() = qual_buf;

                // Adjust read name for splits
                if slice_count > 0 {
                    let mut name_bytes = record_buf.name().map(|n| n.to_vec()).unwrap_or_default();
                    name_bytes.extend_from_slice(format!("_{slice_count}").as_bytes());
                    *record_buf.name_mut() = Some(bstr::BString::from(name_bytes));
                }

                // Update MM/ML tags
                let data = record_buf.data_mut();
                if let Some(mm_val) = new_mm {
                    data.insert(mm_tag.clone(), Value::String(bstr::BString::from(mm_val)));
                } else {
                    data.remove(&mm_tag);
                }
                if let Some(ml_val) = new_ml {
                    data.insert(
                        ml_tag.clone(),
                        Value::Array(noodles::sam::alignment::record_buf::data::field::value::Array::UInt8(ml_val)),
                    );
                } else {
                    data.remove(&ml_tag);
                }

                // Alignment metadata removal (POS, CIGAR, Flags)
                let is_aligned = !record.flags().is_unmapped() || record.alignment_start().is_some();
                if is_aligned {
                    if !warned_about_alignment {
                        println!(
                            "WARNING: Aligned records detected! Removing alignment metadata (CIGAR, mapping position, etc.) to keep BAM records valid after sequence trimming."
                        );
                        warned_about_alignment = true;
                    }
                    *record_buf.flags_mut() = noodles::sam::alignment::record::Flags::UNMAPPED;
                    *record_buf.alignment_start_mut() = None;
                    *record_buf.mapping_quality_mut() = None;
                    *record_buf.cigar_mut() = noodles::sam::alignment::record_buf::Cigar::default();
                    *record_buf.reference_sequence_id_mut() = None;
                    *record_buf.mate_reference_sequence_id_mut() = None;
                    *record_buf.mate_alignment_start_mut() = None;
                    *record_buf.template_length_mut() = 0;
                }

                let group = label_config.create_label(&slice.annotations);

                // Get or create writer for this group
                if !writers.contains_key(&group) {
                    let output_file = format!("{output_folder}/{group}.trimmed.bam");
                    let file = File::create(&output_file).map_err(|err| {
                        let msg = format_output_file_error(&output_file, &err);
                        progress.print_error(msg.clone());
                        progress.clear();
                        anyhow!(msg)
                    })?;
                    let mut writer = noodles::bam::io::Writer::new(file);
                    let output_header = Header::default();
                    writer.write_header(&output_header)?;
                    writers.insert(group.clone(), writer);
                }

                let writer = writers.get_mut(&group).unwrap();
                let output_header = Header::default();
                writer.write_alignment_record(&output_header, &record_buf)?;
                written += 1;
            }

            if written == 0 {
                mark_failed(&progress, &mut failed_trimmed_writer, &read_id);
            } else {
                progress.inc(TRIMMED_IDX);
                if written > 1 {
                    progress.inc(TRIMMED_SPLIT_IDX);
                }
            }
        }

        progress.refresh();
    }

    progress.finish();
    if config.qc.is_active() {
        print_qc_summary_labeled(&config.qc, "Post-trim QC");
    }
    Ok(())
}

