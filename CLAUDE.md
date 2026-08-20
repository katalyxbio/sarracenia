# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Sarracenia is a Rust CLI for pattern-aware demultiplexing of Nanopore/PacBio sequencing reads (FASTQ and BAM, including `MM`/`ML` epigenetic base-modification tags). It finds barcodes/adapters/primers in reads via edit-distance matching (using the `sassy` aligner), classifies the resulting "patterns" of matches, and trims/sorts reads into per-sample output files.

## Build, test, run

```bash
cargo build                          # debug build
RUSTFLAGS="-C target-cpu=native" cargo build --release   # perf build (see .cargo/config.toml)
cargo test                           # run unit tests (many modules have #[cfg(test)] blocks)
cargo test <test_name>               # run a single test
cargo run -- <subcommand> ...        # e.g. cargo run -- kit -k SQK-RBK114-96 -i reads.fastq -o out --maximize
```

- This is a Cargo workspace: root crate (`sarracenia`, lib + bin) plus `benchmarks/` (a separate `publish = false` crate with its own `main.rs` for simulating reads and comparing against other demuxers — not part of the core library).
- `.cargo/config.toml` sets `target-cpu=native` for local builds; `.cargo/config-portable.toml` (swapped in during CI release builds, see `.github/workflows/release.yaml`) targets `x86-64-v3`/`apple-a14` for portable binaries. Don't assume native-only optimizations are safe to rely on in code that ships in releases.
- CI (`release.yaml`) builds with `cargo build --profile dist` and `cargo test --profile dist` (release-like profile with `lto = true`), triggered on `v*` tags.

## Architecture

The pipeline is four stages, each a module under `src/`, invoked as CLI subcommands wired up in `bin/main.rs` (a thin clap `Cli`/`Commands` enum that builds a `*Config` struct from `src/config.rs` and calls into the library):

1. **`annotate`** (`src/annotate/`) — scans reads against query sequences (custom FASTA or a built-in kit preset) and writes a TSV of every match found (`annotator.rs` drives parallel FASTQ/BAM reading via `seq_io`/`noodles` + `rayon`; `searcher.rs`'s `Demuxer` wraps `sassy::Searcher` to find and score barcode/flank candidates; `barcodes.rs` defines query groups and `BarcodeType` (Ftag/Rtag); `cigar_parse.rs`/`interval.rs`/`edit_model.rs` handle CIGAR-derived coordinates, overlap collapsing, and edit-distance cutoffs).
2. **`inspect`** (`src/inspect/`) — summarizes the annotation TSV into "patterns" (see Pattern DSL below), showing how often each pattern occurs across reads.
3. **`filter`** (`src/filter/`) — `pattern.rs` implements the pattern mini-language (parsing pattern strings like `Ftag[fw, *, @left(0..250), >>]` into `PatternElement`/`Cut` structs) and matches them against annotated reads; `filter.rs` reads a user-supplied `filters.txt` of patterns-to-keep and writes only matching rows to the output TSV, populating the `cuts` column that `trim` needs.
4. **`trim`** (`src/trim/`) — `trim.rs` performs the actual read cutting/reverse-complementing/sorting into per-label output files (FASTQ or BAM, optionally gzipped), driven by the `cuts` metadata from `filter`; `mod_recalculator.rs` recomputes BAM `MM`/`ML` tags when a read is trimmed or flipped so base-modification calls stay coordinate-correct.

Orthogonal to those four stages, **`src/qc/`** implements length/mean-quality read filtering. `QcFilter` (thresholds + `evaluate`) is the single source of truth, used from three places: as a **pre-alignment** gate inside `annotate` (both the FASTQ `read_parallel` worker and `process_bam_batch`, reached via `AnnotateConfig.qc` / `KitConfig.qc`), as a **post-trim** gate inside `trim` (both `trim_matches` and `trim_bam_matches`, reached via `TrimConfig.qc` / `KitConfig.trim_qc`), and as the standalone `sarracenia qc` subcommand (`run_qc`, reads in → reads out, no alignment). The two inline gates are separate `QcFilter`s driven by separate CLI flags (`--min-length` etc. vs `--min-trimmed-length` etc.) — a `kit` run can have both active. Read quality is the mean *error probability* converted back to Phred (`-10*log10(mean(10^(-q/10)))`) — the Dorado/NanoFilt/chopper convention — not the arithmetic mean of Phred scores; `mean_error_quality`/`passes_phred` take raw Phred (BAM) and `mean_error_quality_ascii`/`passes_ascii` take a FASTQ `+33` line, so don't mix them up when adding a call site. Drop tallies live in module-level atomics surfaced by `qc_counts()`/`print_qc_summary_labeled()`; because they are process-wide, `annotate` and `trim` each reset them and print their own summary, and any test asserting on them must hold `lock_qc_counts_for_test()`.

Supporting modules:
- **`src/kits/`** — built-in presets. `kits.rs` embeds barcode/adapter sequences for supported Nanopore/PacBio kits at compile time; `use_kit.rs`'s `demux_using_kit` runs the full annotate→filter→trim pipeline for a named kit (the `kit` subcommand, and what `--maximize` affects at the filter-pattern-selection step); `auto_detect.rs` samples reads and picks the best-matching kit automatically.
- **`src/io/`** — shared FASTQ/FASTQ.gz/BAM opening helpers.
- **`src/progress/`** — the whole CLI UI, which is deliberately plain: no ASCII art, no color, no dependency on `indicatif`/`colored`. Informational lines go to stdout via `info()`; each stage's single progress line goes to stderr via `ProgressTracker`, rewritten in place with `\r` on a terminal and appended every 10s when stderr is not a TTY (`Renderer::tty`). The global `--quiet` flag sets a `QUIET` atomic; `info()` and every draw become no-ops, but errors still print. Add a stage by adding a `StageSpec` (name, unit, metric labels) — `metrics[0]` is the running total, and the same labels feed the `--verbose` log file through `log_label()`, so renaming one changes that file's format.

### The Pattern DSL

Central to `inspect`/`filter`/`trim`: a match "pattern" string like:

```
Ftag[fw, *, @left(0..250), >>]__Rtag[<<, rc, *, @right(0..250)]
```

encodes tag type (`Ftag`/`Fflank`/`Rtag`/`Rflank`), orientation (`fw`/`rc`), label matcher (`*`, exact label, or `~substring`), relative position (`@left`/`@right`/`@prev_left` with a distance range), and optional cut markers (`>>`/`<<`, with optional numeric group ids for multi-cut/concat reads). Multiple elements combine with `__`. See the README's "Patterns" section for the full spec — `src/filter/pattern.rs` is the parser/matcher implementation to consult when changing this DSL.

## Notes for changes in this repo

- Kit presets in `src/kits/kits.rs` are large embedded data tables — when adding/fixing a kit, check `data/supported_kits.txt` and the README's kit list for consistency. `test_all_kit_pattern_sets_parse` guards the `LazyLock` pattern sets, which would otherwise only blow up at runtime for whoever ran that kit.
- A kit's `KitConfig.label_config` drives output-file naming through `demux_using_kit` — don't hardcode label options in `use_kit.rs`. Most kits carry the same barcode at both ends and use `only_side: Left`; combinatorial kits (`PacBio-M13`, where the sample is the forward/reverse primer *pair*) use `DOUBLE_LABEL_CONFIG_KEEP_DOUBLE` so both labels reach the filename.
- Pattern sets are per-kit-family, not one-size-fits-all: `SINGLE_LABEL_*` only ever allows `fw` at the right end, `DOUBLE_LABEL_*` is `Ftag`-only, and PacBio needs neither — `PACBIO_96_PATTERNS_*` covers the reverse-complemented trailing barcode of symmetric SMRTbell adapters and `PACBIO_M13_PATTERNS_*` is the only set that references `Rtag`/`Rflank`. Adding a pattern to a shared set changes every kit using it, so prefer a new set.
- Trimming/flipping BAM reads must go through `mod_recalculator.rs` for `MM`/`ML` correctness, and clears alignment metadata (POS/CIGAR/flags/MAPQ) on aligned input — don't bypass this when touching `trim.rs`.
- Paper/reproducibility evaluation code lives in a separate repo (`sarracenia-evals`), not here — `benchmarks/` in this repo is for simulation/comparison tooling only.
