use sarracenia::annotate::annotator::*;
use sarracenia::annotate::barcodes::BarcodeType;
use sarracenia::config::{AnnotateConfig, FilterConfig, KitConfig, QcConfig, TrimConfig};
use sarracenia::filter::filter::filter_from_text_file;
use sarracenia::inspect::inspect;
use sarracenia::kits::use_kit::demux_using_kit;
use sarracenia::progress::progress::{info, set_quiet};
use sarracenia::qc::qc::{QcFilter, run_qc};
use sarracenia::trim::trim::{LabelSide, trim_matches};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Only report errors; suppress progress and informational output
    #[arg(long, global = true, default_value_t = false)]
    quiet: bool,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum UtilsCommands {
    /// Group raw reads by label without trimming
    Pull {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Read FASTQ file (or FASTQ.gz; slower due to unzipping)
        #[arg(short = 'r', long)]
        reads: String,

        /// Output folder path for grouped reads
        #[arg(short = 'o', long)]
        output: String,

        /// Disable label in output filenames
        #[arg(long, default_value_t = false)]
        no_label: bool,

        /// Disable orientation in output filenames
        #[arg(long, default_value_t = false)]
        no_orientation: bool,

        /// Disable flank in output filenames
        #[arg(long, default_value_t = false)]
        no_flanks: bool,

        /// Sort barcode labels in output filenames
        #[arg(long, default_value_t = false)]
        sort_labels: bool,

        /// Only keep left or right label in output filenames
        #[arg(long, conflicts_with = "sort_labels")]
        only_side: Option<LabelSide>,

        /// Write ids of reads with no annotations to this file
        #[arg(long)]
        failed_out: Option<String>,
    },
}

#[derive(Subcommand)]
enum Commands {
    /// Annotate FASTQ/BAM files with barcode information
    Annotate {
        /// Read FASTQ/BAM file (or FASTQ.gz)
        #[arg(short = 'i', long)]
        input: String,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 10)]
        threads: usize,

        /// Output file path
        #[arg(short = 'o', long, default_value = "output.tsv")]
        output: String,

        /// Query files (comma-separated paths)
        #[arg(short = 'q', long, required_unless_present = "kit")]
        queries: Option<String>,

        /// Barcode types (comma-separated: Ftag,Rtag) matching your query file (-q)
        #[arg(short = 'b', long, default_value = "Ftag")]
        barcode_types: String,

        /// Kit name (e.g. SQK-RBK114-24). Conflicts with --queries/--barcode-types
        #[arg(long, conflicts_with = "queries", conflicts_with = "barcode_types")]
        kit: Option<String>,

        /// Flank maximum erors in flank, ONLY set manually when you know what you are doing
        #[arg(long = "flank-max-errors", value_name = "INT")]
        flank_max_errors: Option<usize>,

        /// Enable verbose output for debugging
        #[arg(long, default_value_t = false)]
        verbose: bool,

        /// Barcode: fraction compared to 'perfect' match score for top candidate
        #[arg(long = "min-score", default_value_t = 0.2)]
        min_score: f64,

        /// Barcode: fraction difference between top 2 candidates
        #[arg(long = "min-score-diff", default_value_t = 0.1)]
        min_score_diff: f64,

        /// Also use extended templates (if using kit), i.e. detect fusions, breaks, etc. (slower)
        #[arg(long, default_value_t = false)]
        use_extended: bool,

        /// Edit cost beyond read boundaries
        #[arg(long = "alpha", default_value_t = 0.4)]
        alpha: f32,

        /// Drop reads shorter than this many bases (applied before alignment)
        #[arg(long = "min-length", value_name = "INT")]
        min_length: Option<usize>,

        /// Drop reads longer than this many bases (applied before alignment)
        #[arg(long = "max-length", value_name = "INT")]
        max_length: Option<usize>,

        /// Drop reads with mean read quality below this Phred value (applied before alignment)
        #[arg(long = "min-quality", value_name = "FLOAT")]
        min_quality: Option<f64>,
    },
    /// Filter annotation files based on pattern
    Filter {
        /// Input annotation file
        #[arg(short = 'i', long, required = true)]
        input: String,

        /// Output filtered file path
        #[arg(short = 'o', long, required = true)]
        output: String,

        /// File containing patterns to filter by
        #[arg(short = 'f', long, required = true)]
        file: String,

        /// Write dropped read annotation to this file
        #[arg(long)]
        dropped: Option<String>,

        /// Writes log file (total, kept, dropped)
        #[arg(long, default_value_t = false)]
        verbose: bool,
    },

    /// Trim and sort reads based on filtered annotations
    Trim {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Read FASTQ file (or FASTQ.gz; slower due to unzipping)
        #[arg(short = 'r', long)]
        reads: String,

        /// Output folder path for trimmed reads
        #[arg(short = 'o', long)]
        output: String,

        /// Disable label in output filenames
        #[arg(long, default_value_t = false)]
        no_label: bool,

        /// Disable orientation in output filenames
        #[arg(long, default_value_t = false)]
        no_orientation: bool,

        /// Disable flank in output filenames
        #[arg(long, default_value_t = false)]
        no_flanks: bool,

        /// Sort barcode labels in output filenames
        #[arg(long, default_value_t = false)]
        sort_labels: bool,

        /// Only keep left or right label in output filenames
        #[arg(long, conflicts_with = "sort_labels")]
        only_side: Option<LabelSide>,

        /// Write ids of failed trimmed reads to this file
        #[arg(long)]
        failed_out: Option<String>,

        /// Don't trim reads
        #[arg(long, default_value_t = false)]
        skip_trim: bool,

        /// EXPERIMENTAL, if any Ftag matches rc, reverse complement entire read
        #[arg(long, default_value_t = false)]
        flip: bool,

        /// Writes log file (total, kept, dropped)
        #[arg(long, default_value_t = false)]
        verbose: bool,

        /// Write output FASTQ files as gzip-compressed (.fastq.gz)
        #[arg(long, default_value_t = false)]
        gzip: bool,

        /// Drop trimmed reads shorter than this many bases (applied after trimming)
        #[arg(long = "min-trimmed-length", value_name = "INT")]
        min_trimmed_length: Option<usize>,

        /// Drop trimmed reads longer than this many bases (applied after trimming)
        #[arg(long = "max-trimmed-length", value_name = "INT")]
        max_trimmed_length: Option<usize>,

        /// Drop trimmed reads with mean quality below this Phred value (applied after trimming)
        #[arg(long = "min-trimmed-quality", value_name = "FLOAT")]
        min_trimmed_quality: Option<f64>,
    },

    /// View most common patterns in annotation
    Inspect {
        /// Input filtered annotation file
        #[arg(short = 'i', long)]
        input: String,

        /// Top N
        #[arg(short = 'n', long, default_value_t = 10)]
        top_n: usize,

        /// Write pattern for each read to this file (optional)
        #[arg(short = 'o', long)]
        read_pattern_out: Option<String>,

        /// To summarize results we uses "buckets", such that matches 100 and 103 from the start end up in the same bucket
        #[arg(short = 's', long = "bucket-size", default_value_t = 250)]
        bucket_size: usize,
    },

    /// Run a preset
    Kit {
        /// Kit to use (e.g. SQK-RBK114-24)
        #[arg(short = 'k', long, required_unless_present = "auto_detect")]
        kit: Option<String>,

        /// Auto-detect the barcoding kit from input reads
        #[arg(long, default_value_t = false)]
        auto_detect: bool,

        /// Minimum percentage of sampled reads (0.0 to 100.0) that the best candidate kit must match during auto-detection
        #[arg(long = "auto-detect-threshold", default_value_t = 5.0)]
        auto_detect_threshold: f64,

        /// Input FASTQ file (or FASTQ.gz; slower due to unzipping)
        #[arg(short = 'i', long)]
        input: String,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 10)]
        threads: usize,

        /// Output folder
        #[arg(short = 'o', long)]
        output: String,

        /// Add more 'risky' patterns to demuxing to maximize assigned reads
        #[arg(long, default_value_t = false)]
        maximize: bool,

        /// Enable verbose output for debugging
        #[arg(long, default_value_t = false)]
        verbose: bool,

        /// Fraction compared to 'perfect' match score for top candidate
        #[arg(long = "min-score", default_value_t = 0.2)]
        min_score: f64,

        /// Fraction difference between top 2 candidates
        #[arg(long = "min-score-diff", default_value_t = 0.1)]
        min_score_diff: f64,

        /// Flank maximum errors in flank, ONLY set manually when you know what you are doing
        #[arg(long = "flank-max-errors", value_name = "INT")]
        flank_max_errors: Option<usize>,

        /// Write ids of failed trimmed reads to this file
        #[arg(long)]
        failed_out: Option<String>,

        /// Also use extended templates (if using kit), i.e. detect fusions, breaks, etc. (slower)
        #[arg(long, default_value_t = false)]
        use_extended: bool,

        /// Edit cost beyond read boundaries
        #[arg(long = "alpha", default_value_t = 0.4)]
        alpha: f32,

        /// Write output FASTQ files as gzip-compressed (.fastq.gz)
        #[arg(long, default_value_t = false)]
        gzip: bool,

        /// Drop reads shorter than this many bases (applied before alignment)
        #[arg(long = "min-length", value_name = "INT")]
        min_length: Option<usize>,

        /// Drop reads longer than this many bases (applied before alignment)
        #[arg(long = "max-length", value_name = "INT")]
        max_length: Option<usize>,

        /// Drop reads with mean read quality below this Phred value (applied before alignment)
        #[arg(long = "min-quality", value_name = "FLOAT")]
        min_quality: Option<f64>,

        /// Drop trimmed reads shorter than this many bases (applied after trimming)
        #[arg(long = "min-trimmed-length", value_name = "INT")]
        min_trimmed_length: Option<usize>,

        /// Drop trimmed reads longer than this many bases (applied after trimming)
        #[arg(long = "max-trimmed-length", value_name = "INT")]
        max_trimmed_length: Option<usize>,

        /// Drop trimmed reads with mean quality below this Phred value (applied after trimming)
        #[arg(long = "min-trimmed-quality", value_name = "FLOAT")]
        min_trimmed_quality: Option<f64>,
    },

    /// Filter reads by length and mean quality, without any alignment
    Qc {
        /// Input FASTQ/BAM file (or FASTQ.gz)
        #[arg(short = 'i', long)]
        input: String,

        /// Output file for reads that pass (same format as the input)
        #[arg(short = 'o', long)]
        output: String,

        /// Also write reads that fail to this file
        #[arg(long)]
        dropped: Option<String>,

        /// Drop reads shorter than this many bases
        #[arg(long = "min-length", value_name = "INT")]
        min_length: Option<usize>,

        /// Drop reads longer than this many bases
        #[arg(long = "max-length", value_name = "INT")]
        max_length: Option<usize>,

        /// Drop reads with mean read quality below this Phred value.
        /// Quality is -10*log10(mean error probability), as used by Dorado/NanoFilt/chopper
        #[arg(long = "min-quality", value_name = "FLOAT")]
        min_quality: Option<f64>,

        /// Write output FASTQ files as gzip-compressed (.fastq.gz)
        #[arg(long, default_value_t = false)]
        gzip: bool,

        /// Writes log file (total, kept, dropped)
        #[arg(long, default_value_t = false)]
        verbose: bool,
    },
}

fn main() {
    // Make sure that AVX2 is supported by the current CPU if it was used during compilation.
    ensure_simd::ensure_simd();

    let cli = Cli::parse();

    set_quiet(cli.quiet);

    match &cli.command {
        Commands::Annotate {
            input,
            threads,
            output,
            queries,
            barcode_types,
            kit,
            flank_max_errors,
            verbose,
            min_score,
            min_score_diff,
            use_extended,
            alpha,
            min_length,
            max_length,
            min_quality,
        } => {
            info("Starting annotation...");
            let annotate_config = AnnotateConfig {
                max_flank_errors: *flank_max_errors,
                alpha: *alpha,
                n_threads: *threads as u32,
                verbose: *verbose,
                min_score: *min_score,
                min_score_diff: *min_score_diff,
                use_extended: *use_extended,
                qc: QcFilter::new(*min_length, *max_length, *min_quality),
            };

            if let Some(kit_name) = kit.as_ref() {
                match annotate_with_kit(input, output, kit_name.as_str(), &annotate_config) {
                    Ok(_) => info("Annotation complete!"),
                    Err(e) => println!("Error during processing: {e}"),
                }
                return;
            }

            // Split comma-separated query paths into Vec<String>
            let queries_value = queries
                .as_ref()
                .expect("--queries is required unless --kit is provided");
            let query_files: Vec<String> = queries_value
                .split(',')
                .map(|s| s.trim().to_string())
                .collect();

            let query_files_refs: Vec<&str> = query_files.iter().map(|s| s.as_str()).collect();

            // Parse barcode types
            let barcode_types_vec: Vec<BarcodeType> = barcode_types
                .split(',')
                .map(|s| match s.trim() {
                    "Ftag" => BarcodeType::Ftag,
                    "Rtag" => BarcodeType::Rtag,
                    _ => {
                        panic!("Unknown barcode type: {s}, use one of: Ftag, Rtag")
                    }
                })
                .collect();

            match annotate_with_files(
                input,
                query_files_refs,
                barcode_types_vec,
                output,
                &annotate_config,
            ) {
                // Convert fractions to raw scores
                Ok(_) => info("Annotation complete!"),
                Err(e) => println!("Error during processing: {e}"),
            }
        }

        Commands::Filter {
            input,
            output,
            file,
            dropped,
            verbose,
        } => {
            info("Starting filtering...");
            let filter_config = FilterConfig { verbose: *verbose };

            match filter_from_text_file(input, file, output, dropped.as_deref(), &filter_config) {
                Ok(_) => info("Filtering successful!"),
                Err(e) => println!("Filtering failed: {e}"),
            }
        }

        Commands::Trim {
            input,
            reads,
            output,
            no_label,
            no_orientation,
            no_flanks,
            sort_labels,
            only_side,
            failed_out,
            skip_trim,
            flip,
            verbose,
            gzip,
            min_trimmed_length,
            max_trimmed_length,
            min_trimmed_quality,
        } => {
            info("Starting trimming...");
            let trim_config = TrimConfig {
                add_labels: !no_label,
                add_orientation: !no_orientation,
                add_flank: !no_flanks,
                sort_labels: *sort_labels,
                only_side: *only_side,
                failed_trimmed_writer: failed_out.clone(),
                write_full_header: true, // Maybe make this optional but dont see a reason why you would not want this
                skip_trim: *skip_trim,
                flip: *flip,
                verbose: *verbose,
                gzip: *gzip,
                qc: QcFilter::new(
                    *min_trimmed_length,
                    *max_trimmed_length,
                    *min_trimmed_quality,
                ),
            };
            match trim_matches(input, reads, output, &trim_config) {
                Ok(_) => info("Trimming complete!"),
                Err(e) => println!("Trimming failed: {e}"),
            }
        }

        Commands::Inspect {
            input,
            top_n,
            read_pattern_out,
            bucket_size,
        } => {
            info("Inspecting...");

            match inspect::inspect(input, *top_n, read_pattern_out.clone(), *bucket_size) {
                Ok(_) => info("Inspection complete!"),
                Err(e) => println!("Inspection failed: {e}"),
            }
        }

        Commands::Kit {
            kit,
            auto_detect,
            auto_detect_threshold,
            input,
            threads,
            output,
            maximize,
            verbose,
            min_score,
            min_score_diff,
            flank_max_errors,
            failed_out,
            use_extended,
            alpha,
            gzip,
            min_length,
            max_length,
            min_quality,
            min_trimmed_length,
            max_trimmed_length,
            min_trimmed_quality,
        } => {
            let kit_name = if *auto_detect {
                match sarracenia::kits::auto_detect::auto_detect_kit(input.as_str(), *auto_detect_threshold) {
                    Ok(detected) => detected,
                    Err(e) => {
                        println!("Auto-detection failed: {e}");
                        return;
                    }
                }
            } else {
                kit.clone().expect("Kit name is required unless auto-detect is enabled")
            };

            let kit_config = KitConfig {
                kit_name,
                threads: *threads,
                output_folder: output.clone(),
                maximize: *maximize,
                verbose: *verbose,
                min_score: *min_score,
                min_score_diff: *min_score_diff,
                max_flank_errors: *flank_max_errors,
                failed_out: failed_out.clone(),
                use_extended: *use_extended,
                alpha: *alpha,
                gzip: *gzip,
                qc: QcFilter::new(*min_length, *max_length, *min_quality),
                trim_qc: QcFilter::new(
                    *min_trimmed_length,
                    *max_trimmed_length,
                    *min_trimmed_quality,
                ),
            };

            if let Err(e) = demux_using_kit(input.as_str(), &kit_config) {
                println!("Demultiplexing failed: {e}");
            }
        }

        Commands::Qc {
            input,
            output,
            dropped,
            min_length,
            max_length,
            min_quality,
            gzip,
            verbose,
        } => {
            info("Starting read QC...");
            let qc_config = QcConfig {
                filter: QcFilter::new(*min_length, *max_length, *min_quality),
                gzip: *gzip,
                verbose: *verbose,
            };

            match run_qc(input, output, dropped.as_deref(), &qc_config) {
                Ok(_) => info("Read QC complete!"),
                Err(e) => println!("Read QC failed: {e}"),
            }
        }
    }
}
