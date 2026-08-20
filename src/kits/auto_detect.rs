use crate::annotate::barcodes::BarcodeGroup;
use crate::annotate::edit_model::get_edit_cut_off;
use crate::annotate::searcher::Demuxer;
use crate::io::io::open_fastq;
use crate::progress::progress::info;
use anyhow::{anyhow, Result};
use seq_io::fastq::Record;
use std::fs::File;

const CANDIDATE_KITS: &[&str] = &[
    "SQK-RBK114-96",
    "SQK-RBK114-24",
    "SQK-NBD114-96",
    "SQK-NBD114-24",
    "SQK-16S114-24",
    "SQK-LWB001",
    "SQK-PCB114-24",
    "EXP-NBD104",
    "EXP-NBD114",
    "EXP-PBC001",
    "EXP-PBC096",
    "SQK-RAB204",
    "SQK-RBK001",
    "SQK-RBK111-96",
    "SQK-RBK111-24",
    "SQK-RBK004",
    "SQK-RLB001",
    "SQK-RPB114-24",
    "VSK-VMK001",
    "VSK-VMK004",
    "SQK-MAB114-24",
    "PacBio-M13",
    "PacBio-96A",
];

/// Auto-detects the most likely barcoding kit by matching the first 5000 reads against
/// representative kits, weighted by the kit's flanking sequence length.
pub fn auto_detect_kit(read_file: &str, min_match_pct: f64) -> Result<String> {
    info("Auto-detecting barcoding kit from input file...");

    // 1. Read first 5000 reads
    let mut sample_reads = Vec::new();
    if read_file.to_ascii_lowercase().ends_with(".bam") {
        let mut reader = File::open(read_file)
            .map(noodles::bam::io::Reader::new)
            .map_err(|e| anyhow!("Failed to open BAM file '{read_file}': {e}"))?;
        let _header = reader.read_header()?;
        for result in reader.records() {
            let record = result.map_err(|e| anyhow!("Failed to read BAM record: {e}"))?;
            let seq = record.sequence();
            let mut seq_bytes = Vec::with_capacity(seq.len());
            for base in seq.iter() {
                seq_bytes.push(base.to_ascii_uppercase());
            }
            sample_reads.push(seq_bytes);
            if sample_reads.len() >= 5000 {
                break;
            }
        }
    } else {
        let mut reader = open_fastq(read_file);
        while let Some(record) = reader.next() {
            let seqrec = record.map_err(|e| anyhow!("Failed to read FASTQ record: {e}"))?;
            sample_reads.push(seqrec.seq().to_vec());
            if sample_reads.len() >= 5000 {
                break;
            }
        }
    }

    if sample_reads.is_empty() {
        return Err(anyhow!("Input file has no reads, cannot auto-detect kit."));
    }

    info(format!(
        "Sampled {} reads for auto-detection.",
        sample_reads.len()
    ));

    // 2. Score candidate kits
    let mut best_kit = None;
    let mut max_score = 0;
    let mut best_match_pct = 0.0;

    for &kit_name in CANDIDATE_KITS {
        // Load the barcode groups for the candidate kit
        // PacBio presets are supported too
        let mut total_flank_len = 0;
        let mut groups = BarcodeGroup::new_from_kit(kit_name, false);
        for group in &mut groups {
            total_flank_len += group.flank_prefix.len() + group.flank_suffix.len();
            let edit_cut_off = get_edit_cut_off(group.get_effective_len());
            group.set_flank_threshold(edit_cut_off);
        }

        let mut demuxer = Demuxer::new(0.4, false, 0.2, 0.1);
        for group in groups {
            demuxer.add_query_group(group);
        }

        let mut match_count = 0;
        for read in &sample_reads {
            let matches = demuxer.demux("temp", read);
            if !matches.is_empty() {
                match_count += 1;
            }
        }

        let score = match_count * std::cmp::max(total_flank_len, 1);
        let match_pct = (match_count as f64 / sample_reads.len() as f64) * 100.0;

        if match_count > 0 {
            info(format!(
                "  Kit {} matched {}/{} reads ({:.2}%, flank_len={}, score={}).",
                kit_name,
                match_count,
                sample_reads.len(),
                match_pct,
                total_flank_len,
                score
            ));
        }

        if score > max_score {
            max_score = score;
            best_kit = Some(kit_name.to_string());
            best_match_pct = match_pct;
        }
    }

    if let Some(kit) = best_kit {
        if best_match_pct >= min_match_pct {
            info(format!(
                "Auto-detected kit: {} with score {} (match rate: {:.2}%)",
                kit, max_score, best_match_pct
            ));
            Ok(kit)
        } else {
            Err(anyhow!(
                "Could not confidently auto-detect the kit. Best candidate '{}' only matched {:.2}% of reads, which is below the threshold of {:.2}%. Please specify the kit manually using -k.",
                kit, best_match_pct, min_match_pct
            ))
        }
    } else {
        Err(anyhow!(
            "No candidate barcoding kit matched any of the sampled reads. Threshold is {:.2}%. Please specify the kit manually using -k.",
            min_match_pct
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kits::kits::lookup_barcode_seq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_auto_detect_rbk114_96() {
        // Construct a realistic read seq for SQK-RBK114-96 using RBK01
        let front = "GCTTGGGTGTTTAACC";
        let rear = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
        let barcode_seq = lookup_barcode_seq("RBK01").expect("Failed to get RBK01 sequence");
        let full_read = format!("{front}{barcode_seq}{rear}ACGTACGTACGTACGTACGT");

        let mut temp_file = NamedTempFile::new().unwrap();
        for i in 0..10 {
            writeln!(
                temp_file,
                "@read_{}\n{}\n+\n{}",
                i,
                full_read,
                "F".repeat(full_read.len())
            )
            .unwrap();
        }

        let detected = auto_detect_kit(temp_file.path().to_str().unwrap(), 5.0).unwrap();
        assert_eq!(detected, "SQK-RBK114-96");
    }
}
