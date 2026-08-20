use crate::qc::qc::QcFilter;
use crate::trim::trim::LabelSide;

#[derive(Debug, Clone)]
pub struct AnnotateConfig {
    pub max_flank_errors: Option<usize>,
    pub alpha: f32,
    pub n_threads: u32,
    pub verbose: bool,
    pub min_score: f64,
    pub min_score_diff: f64,
    pub use_extended: bool,
    /// Length/quality gate applied to raw reads before any alignment work.
    pub qc: QcFilter,
}

/// Standalone `sarracenia qc`: reads in, surviving reads out, no alignment.
#[derive(Debug, Clone)]
pub struct QcConfig {
    pub filter: QcFilter,
    pub gzip: bool,
    pub verbose: bool,
}

#[derive(Debug, Clone)]
pub struct FilterConfig {
    pub verbose: bool,
}

#[derive(Debug, Clone)]
pub struct TrimConfig {
    pub add_labels: bool,
    pub add_orientation: bool,
    pub add_flank: bool,
    pub sort_labels: bool,
    pub only_side: Option<LabelSide>,
    pub failed_trimmed_writer: Option<String>,
    pub write_full_header: bool,
    pub skip_trim: bool,
    pub flip: bool,
    pub verbose: bool,
    pub gzip: bool,
    /// Length/quality gate applied to each fragment *after* trimming, just before it
    /// is written out. Independent of `AnnotateConfig.qc`, which sees the raw read.
    pub qc: QcFilter,
}

#[derive(Debug, Clone)]
pub struct KitConfig {
    pub kit_name: String,
    pub threads: usize,
    pub output_folder: String,
    pub maximize: bool,
    pub verbose: bool,
    pub min_score: f64,
    pub min_score_diff: f64,
    pub max_flank_errors: Option<usize>,
    pub failed_out: Option<String>,
    pub use_extended: bool,
    pub alpha: f32,
    pub gzip: bool,
    /// Length/quality gate applied to raw reads before any alignment work.
    pub qc: QcFilter,
    /// Length/quality gate applied to trimmed reads just before they are written.
    pub trim_qc: QcFilter,
}

