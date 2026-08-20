use crate::annotate::annotator::annotate_with_kit;
use crate::config::{AnnotateConfig, FilterConfig, KitConfig, TrimConfig};
use crate::filter::filter::filter;
use crate::inspect::inspect::inspect;
use crate::kits::kits::*;
use crate::trim::trim::trim_matches;
use anyhow::anyhow;
use crate::progress::progress::{info, is_quiet};
use std::path::Path;

pub fn demux_using_kit(fastq_file: &str, config: &KitConfig) -> anyhow::Result<()> {
    let kit_name = config.kit_name.as_str();
    let output_folder = config.output_folder.as_str();
    // Create output folder if not exists yet
    if !Path::new(output_folder).exists() {
        std::fs::create_dir_all(output_folder)?;
    }

    let kit_info = get_kit_info(kit_name);

    // Print some kit info
    if !is_quiet() {
        info("");
        info("Kit info");
        info(format!("Kit name: {}", kit_info.name));
        info(format!(
            "Kit type: {}",
            if config.maximize { "Maximize" } else { "Safe" }
        ));
        for tmpl in kit_info.templates {
            info(format!(
                "Barcodes: {} - {}",
                tmpl.barcodes.from, tmpl.barcodes.to
            ));
        }
        if config.qc.is_active() {
            info(format!("Read QC (pre-alignment): {}", config.qc.describe()));
        }
        if config.trim_qc.is_active() {
            info(format!("Read QC (post-trim): {}", config.trim_qc.describe()));
        }
    }

    info("");
    info("Annotating reads...");
    let annotate_config = AnnotateConfig {
        max_flank_errors: config.max_flank_errors,
        alpha: config.alpha,
        n_threads: config.threads as u32,
        verbose: config.verbose,
        min_score: config.min_score,
        min_score_diff: config.min_score_diff,
        use_extended: config.use_extended,
        qc: config.qc,
    };
    annotate_with_kit(
        fastq_file,
        format!("{output_folder}/annotation.tsv").as_str(),
        kit_name,
        &annotate_config,
    )?;

    // // After annotating we show inspect
    if !is_quiet() {
        info("");
        info("Top 10 most common patterns");
        let pattern_per_read_out = format!("{output_folder}/pattern_per_read.tsv");
        inspect(
            format!("{output_folder}/annotation.tsv").as_str(),
            10,
            Some(pattern_per_read_out),
            250,
        )
        .map_err(|e| anyhow!("{e}"))?;
        info(format!(
            "Want to see more patterns? Run: `sarracenia inspect {output_folder}/annotation.tsv -n 100`"
        ));
    }

    // Filter
    info("");
    info("Filtering reads...");

    let patterns = if config.maximize {
        (kit_info.maximize_patterns)()
    } else {
        (kit_info.safe_patterns)()
    };
    let filter_config = FilterConfig {
        verbose: config.verbose,
    };

    filter(
        format!("{output_folder}/annotation.tsv").as_str(),
        format!("{output_folder}/filtered.tsv").as_str(),
        None,
        patterns,
        &filter_config,
    )
    .map_err(|e| anyhow!("{e}"))?;

    // Trimming
    info("");
    info("Trimming reads...");
    // Naming comes from the kit itself: most kits carry the same barcode at both ends
    // and want just the left label, but a combinatorial kit needs both.
    let label_config = kit_info.label_config;
    let trim_config = TrimConfig {
        add_labels: label_config.include_label,
        add_orientation: label_config.include_orientation,
        add_flank: label_config.include_flank,
        sort_labels: label_config.sort_labels,
        only_side: label_config.only_side,
        failed_trimmed_writer: config.failed_out.clone(),
        write_full_header: true,
        skip_trim: false,
        flip: false,
        verbose: config.verbose,
        gzip: config.gzip,
        qc: config.trim_qc,
    };
    trim_matches(
        format!("{output_folder}/filtered.tsv").as_str(),
        fastq_file,
        output_folder,
        &trim_config,
    )?;

    info("");
    info("Done!");
    Ok(())
}
