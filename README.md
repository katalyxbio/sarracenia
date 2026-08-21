# Sarracenia — Pattern aware demux

![Sarracenia](resources/sarracenia.png)

## Why Sarracenia?

- **>1000× fewer trimming errors** compared to Dorado.
- Equivalent or **better assemblies**.
- **Contamination-free** assemblies by removing artefact reads.
- **High-performance Rust-native BAM support** via Noodles.
- **Dynamic Epigenetic Base Modification Recalculation** (`MM`/`ML` tags) when reads are trimmed or flipped.
- **Embedded PacBio presets** (SMRTbell & M13) with zero external database dependencies.
- Easily applicable to **custom experiments**.
- Still **very fast**.

If you have any issues or if something is unclear, just create an [issue](https://github.com/katalyxbio/sarracenia/issues).



## Quick links

- [Installing Sarracenia](#installing-sarracenia)
- [Quickstart](#quickstart)
- [Read QC: length & quality filtering](#read-qc-length--quality-filtering)
- [BAM & Epigenetic Base Modification Handling](#bam--epigenetic-base-modification-handling)
- [In-depth inspection of Nanopore & PacBio kit results](#in-depth-inspection-of-nanopore--pacbio-kit-results)
  - [Annotate](#annotate)
  - [Inspect](#inspect)
  - [Filter](#filter)
  - [Trim](#trim)
  - [The `--maximize` flag explained in more detail](#themaximize-flag-explained-in-more-detail-kit-command-only)
- [Custom experiment](#custom-experiment)
  - [Creating a query Fasta](#creating-a-query-fasta)
  - [Single end](#single-end)
  - [Dual end](#dual-end)
- [Custom experiment with mixed sequences](#custom-experiment-with-mixed-sequences)
- [Output columns (annotate & filter)](#output-columns-annotate--filter)
- [Patterns](#patterns)
  - [How to handle concat reads](#how-to-handle-concat-reads)
- [Paper evals](#paper-evals)
- [Notes & tips](#notes--tips)
- [License](#license)


## Installing Sarracenia
Sarracenia is written in Rust.

### Executables
Download the latest release for your platform from
[releases](https://github.com/katalyxbio/sarracenia/releases), then make it executable and
put it somewhere on your `PATH`:

``` bash
chmod +x sarracenia-x86_64-unknown-linux-gnu
mv sarracenia-x86_64-unknown-linux-gnu ~/.local/bin/sarracenia
```

Each release also ships conda packages for the same three platforms, if you would rather
install into an environment:

``` bash
conda install ./sarracenia-0.3.3-linux-64.conda
```

### From source (recommended)

> [!IMPORTANT]
> Current git version is around 5x faster than the releases by using [Sassy V2](https://github.com/RagnarGrootKoerkamp/sassy) and other speedups. But 
> results can slightly differ (<0.01% or so).

Check whether Rust is installed:

```bash
rustc --version
```

If not install it via [rustup](https://rustup.rs/), more info in their
[docs](https://rust-lang.github.io/rustup/installation/index.html). Use `rustup
update` to get the latest stable version.

Then clone the repository and build it:

```bash
git clone https://github.com/katalyxbio/sarracenia
cd sarracenia
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

The binary lands in `target/release/sarracenia`; copy it onto your `PATH` to use it from
anywhere.

See [here](https://github.com/ragnargrootkoerkamp/ensure_simd) for
details on `target-cpu=native`.



## Quickstart

Sarracenia includes built-in kit *presets* for many Nanopore and PacBio kits. Presets let you run analyses quickly, but we recommend reading **Understanding Sarracenia** to interpret the results correctly.

Basic command:

```bash
sarracenia kit -k <kit-name> -i <reads.fastq|reads.bam> -o <output_folder> --maximize
```

The `--maximize` option is recommended (e.g., for assembly) unless you need an ultra-strict barcode configuration ([see here](#themaximize-flag-explained-in-more-detail-kit-command-only) for more details). You can also pass `--gzip` to write 
`fastq.gz` files when using FASTQ input, but note that the zipping has a performance penalty.

### Kit Auto-detection
If you do not know which kit was used for your sequencing run, you can let Sarracenia auto-detect it by passing the `--auto-detect` flag instead of specifying `-k`:
```bash
sarracenia kit --auto-detect -i reads.fastq -o output_folder --maximize
```
This will sample the first 100 reads of your input file (works with both FASTQ and BAM), align them against all representative preset kits (including Nanopore and PacBio kits), and automatically choose the kit with the highest number of confident matches to run the demultiplexing pipeline!

### Native barcoding example (SQK-NBD114-96)

```bash
sarracenia kit -k SQK-NBD114-96 -i reads.fastq -o output_folder --maximize
```

This uses a conservative flank-based edit-distance cutoff. If many reads are missed during `annotate`, you can relax the flank error threshold, for example:

```bash
--flank-max-errors 5
```

—but always inspect the results after changing thresholds to avoid random matches (which show up as `Fflank` matches).

### Rapid barcoding example (SQK-RBK114-96)

```bash
sarracenia kit -k SQK-RBK114-96 -i reads.fastq -o output_folder --maximize
```

### PacBio Kits (Sequel II & M13 Plate)
Sarracenia embeds presets for PacBio kits at compile-time. Available PacBio presets include:
- `PacBio-M13`: Combinatorial plate with 16 forward barcodes (`bc1002`–`bc1017`) and 24 reverse barcodes (`bc1050`–`bc1073`), generating 384 unique barcode pairs.
- `PacBio-96A`: Indexed SMRTbell adapters for plate 96A (`bc2001`–`bc2096`).
- `PacBio-96B`: Indexed SMRTbell adapters for plate 96B (`bc2097`–`bc2192`).
- `PacBio-96C`: Indexed SMRTbell adapters for plate 96C (`bc2193`–`bc2288`).
- `PacBio-96D`: Indexed SMRTbell adapters for plate 96D (`bc2289`–`bc2384`).

Usage:
```bash
sarracenia kit -k PacBio-M13 -i reads.bam -o output_folder --maximize
```

The SMRTbell adapters are symmetric, so a full-length `PacBio-96*` read carries the same
barcode at both ends with the trailing one reverse-complemented; reads barcoded on only
one end are accepted too.

`PacBio-M13` is **combinatorial** — the sample is the *pair* of primers, not either one
alone — so its output files are named after both, e.g.
`M13_bc1002_F__M13_bc1050_R.trimmed.bam`. The safe pattern set only keeps reads where
both primers were found. `--maximize` additionally keeps reads where just one was
resolved, under a single-label name (`M13_bc1002_F.trimmed.bam`) so they are easy to tell
apart from properly paired samples — those reads cannot be assigned to one of the 384
pairs, so don't treat such a file as a sample.

For a list of all supported kits, check `data/supported_kits.txt`. 



## Read QC: length & quality filtering

Sarracenia can drop reads that are too short, too long, or too low quality. This works in
three places, using the same thresholds and the same code:

**1. Inline, as a pre-alignment gate.** Add the flags to `kit` or `annotate` and failing
reads are discarded *before* any barcode alignment happens, so you don't pay for aligning
reads you were going to throw away:

```bash
sarracenia kit -k SQK-NBD114-96 -i reads.fastq -o output_folder --maximize \
    --min-length 1000 --min-quality 10
```

**2. Inline, as a post-trim gate.** The `--min-trimmed-*` flags on `kit` or `trim` judge
each fragment *after* its barcodes and adapters have been cut off, so what you filter on
is what actually lands in the output file:

```bash
sarracenia kit -k SQK-NBD114-96 -i reads.fastq -o output_folder --maximize \
    --min-trimmed-length 1000 --min-trimmed-quality 10
```

**3. Standalone, via `sarracenia qc`.** Reads in, surviving reads out, in the same format as
the input (FASTQ, FASTQ.gz or BAM). No alignment is performed:

```bash
sarracenia qc -i reads.fastq -o clean.fastq --min-length 1000 --min-quality 10
```

Flags:

```
  # raw read, before alignment (kit, annotate, qc)
  --min-length <INT>             Drop reads shorter than this many bases
  --max-length <INT>             Drop reads longer than this many bases
  --min-quality <FLOAT>          Drop reads with mean read quality below this Phred value

  # trimmed fragment, before it is written (kit, trim)
  --min-trimmed-length <INT>     Drop trimmed reads shorter than this many bases
  --max-trimmed-length <INT>     Drop trimmed reads longer than this many bases
  --min-trimmed-quality <FLOAT>  Drop trimmed reads with mean quality below this Phred value
```

The two gates are independent and can be combined. A 1100 bp read carrying a 90 bp
adapter+barcode passes `--min-length 1000` but its 1010 bp insert is what
`--min-trimmed-length` sees, and a read whose adapters are high quality can still fail
`--min-trimmed-quality` if the insert itself is poor. A read whose every fragment is
rejected post-trim is counted as failed and, if `--failed-out` is set, its ID is written
there alongside the reads that could not be trimmed at all.

`sarracenia qc` additionally accepts `--dropped <file>` to write the rejected reads, and
`--gzip` to compress FASTQ output (output paths ending in `.gz` are compressed
automatically).

### How read quality is computed

Mean **error probability**, converted back to Phred:

```
Q = -10 * log10( mean( 10^(-q_i/10) ) )
```

This is the same definition used by Dorado, NanoFilt, chopper and filtlong, so
`--min-quality 10` means what you expect coming from those tools. Note it is deliberately
*not* the arithmetic mean of the Phred scores — that overestimates quality whenever a read
carries a few very bad bases. A read of four Q40 bases and one Q0 base scores **Q7**, not
Q32.

The same formula is used for both gates; the post-trim gate simply scores the trimmed
quality string rather than the full one.

> [!NOTE]
> `--min-length`/`--max-length`/`--min-quality` are measured on the **raw** read, before
> trimming. Use the `--min-trimmed-*` flags if you want the thresholds applied to the
> reads that actually end up in the output.

> [!NOTE]
> BAM records with no quality string (`*`) pass the quality check rather than being
> silently dropped; Sarracenia prints a warning reporting how many such reads it saw.



## BAM & Epigenetic Base Modification Handling

Sarracenia natively supports BAM file formats using high-performance, Rust-native `Noodles` integration.

### Epigenetic Base Modifications (MM/ML Tags)
When trimming or reverse-complementing/flipping BAM reads:
- Sarracenia dynamically recalculates `MM` (base modification type/strand) and `ML` (base modification probabilities) tags.
- The recalculation accurately adjusts coordinate-based index skips and mirrors probability values depending on trimmed coordinate segments and strand direction.

### Automatic Alignment Clearance
Trimming/flipping invalidates original genomic alignments:
- If you supply an aligned BAM as input to `sarracenia trim`, Sarracenia will **automatically remove alignment metadata** (POS, CIGAR, flags, MAPQ) to guarantee valid BAM outputs.
- A yellow warning is printed to the console during execution if alignment metadata is cleared.



## In-depth inspection of Nanopore & PacBio kit results

Sarracenia's typical manual workflow: **annotate → inspect → filter → trim**.

### Annotate

Run `annotate` to find matches in reads and output an annotation table (TSV):

```bash
sarracenia annotate --kit SQK-RBK114-96 -i pass_sample.fastq -t 10 -o anno.tsv
```

(Using 10 threads.)

Example `anno.tsv` rows:

```
read_id read_len        rel_dist_to_end read_start_bar  read_end_bar    read_start_flank        read_end_flank   bar_start       bar_end match_type  flank_cost      barcode_cost    label   strand  cuts
c5f925b2-fc0b-4053-b615-d70950d41436    19783   14      14      104     14      104     0       0       Fflank   14      14      flank   Fwd
dbca7fb9-d6c8-4417-8ae7-bc32ebce9b27    2972    29      48      70      29      111     0       23      Ftag     13      7       BC29    Fwd
6c089f0a-50cd-4215-94f0-c7babb87f5fe    7599    27      51      74      27      121     0       23      Ftag     11      5       BC45    Fwd
..etc
```

For column descriptions see [Output columns (annotate & filter)](#output-columns-annotate--filter) below.

This file shows, per read, which barcodes/flanks were matched and with what costs. Use `inspect` next to summarize patterns.

### Inspect

Summarize patterns across the annotation file:

```bash
sarracenia inspect -i anno.tsv
```

By default `inspect` shows the top 10 pattern groups; use `-n <amount>` to increase.

Example summary:

```
Found 64 unique patterns
  Pattern 1: 82421 occurrences
    Ftag[fw, *, @left(0..250)]
  Pattern 2: 5003 occurrences
    Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @right(0..250)]
  Pattern 3: 3545 occurrences
    Fflank[fw, *, @left(0..250)]
  ...
Showed 10 / 64 patterns
Inspection complete!
```

Some observations:
- `Ftag` on the left (`@left`) is the expected pattern for the rapid barcoding kit - that is good news.
- A contamination pattern can be a barcode on the left *and* another barcode on the right (`@right`). We can decide to just trim of the right side (see filtering later)
- `Fflank` indicates that flanks matched but no confident barcode was found.
- `@prev_left` indicates additional tags close to a previous element (e.g., double-barcode ligation).

#### Per-read patterns

To output the selected pattern per read:

```bash
sarracenia inspect -i anno.tsv -o pattern_per_read.tsv
```

Example `pattern_per_read.tsv` contents:

```
85ef... \t Ftag[fw, *, @left(0..250)]
2f67... \t Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250)]
...
```

This is useful when you want to inspect a single "weird" read in detail.

### Filter

Create a `filters.txt` file listing the patterns you want to keep, one per line.
For example:

```
Ftag[fw, *, @left(0..250)]
Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @right(0..250)]
Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250)]
```

Then run:

```bash
sarracenia filter -i anno.tsv -f filters.txt -o filtered.tsv
```

The resulting `filtered.tsv` contains only reads that match the specified patterns.

#### Cutting / trimming metadata

Sarracenia needs to know *where* to cut reads for trimming. You mark cut positions by adding `>>` (cut **after** this element) or `<<` (cut **before** this element) inside the tag's bracket list. Where within the brackets does not matter.

Examples (note the comma-separated fields inside the brackets):

```text
Ftag[fw, *, @left(0..250), >>]
Ftag[fw, *, @left(0..250), >>]__Ftag[<<, fw, *, @right(0..250)]
Ftag[fw, *, @left(0..250)]__Ftag[fw, *, @prev_left(0..250), >>]
```

In the middle pattern we retain the read sequence between the left tag (cut after it) and the right tag (cut before it).

Run `filter` again (same command as above) to populate the `cuts` column in `filtered.tsv`. This is required before trimming.



## Trim

Trim reads using the `cuts` metadata produced by `filter`:

```bash
sarracenia trim -i filtered.tsv -r reads.fastq -o trimmed
```
You can also pass `--gzip` to write 
`fastq.gz` files but note that the zipping has a performance penalty. If the input is a BAM file, Sarracenia automatically outputs matched BAM files (e.g. `BC01.trimmed.bam`).

Output files are organized by pattern-based folder/filenames, for example:

```
BC14_fw__BC14_fw.trimmed.fastq   BC31_fw__BC04_fw.trimmed.fastq  ...
```

If you prefer different filename conventions, use these flags:

```
  --no-label               Disable label in output filenames
  --no-orientation         Disable orientation in output filenames
  --no-flanks              Disable flanks in output filenames
  --sort-labels            Sort barcode labels in output filenames
  --only-side <left|right> Only keep left or right label in output filenames
```

Example to remove orientation and keep only the left label:

```bash
sarracenia trim -i filtered.tsv -r reads.fastq -o trimmed --no-orientation --only-side left
```

Gives:
```
BC01.trimmed.fastq  BC11.trimmed.fastq ...
```

### The `--maximize` flag explained in more detail (`kit` command only)

This flag affects the *filter* step. After locating barcodes/flanks, we use them to assign each read to a sample and perform trimming.  
If multiple barcodes/flanks are detected, we apply _pattern filters_ to decide which reads pass. 
I.e. for which reads we are confident enough to assign a sample despite multiple hits.

In the default (e.g. `safe`) mode, only reads with an unambiguous barcode ligation pattern are allowed to pass.  
For example, in the rapid barcoding kit:

- (1) `[BC1]` (assigned: `BC1`)  
- (2) `[BC1][BC1]` (assigned: `BC1`)  

Here, `(1)` is what we expect from the experimental setup, while `(2)` is also commonly observed (see [paper](https://www.biorxiv.org/content/10.1101/2025.10.22.683865v1)).  
Since both barcodes in `(2)` are the same (BC01), we assign the read to BC01 and trim both.

For applications like assembly, you _might_ want to retain more reads, even at the cost of introducing some errors. This is where `--maximize` comes in.  
In addition to the `safe` patterns above, it also allows (among others):

- (1) `[BC2][BC1]` (assigned: BC1)  
- (2) `[BC1][BC2]` (assigned: BC1)  

In `(1)`, we assume BC1 was ligated first (creating `[BC1]-`), followed by another ligation _after_ pooling, resulting in `[BC2][BC1]`.  
The original sample is therefore likely BC1, and we assign it accordingly.  

In `(2)`, we observe an unexpected barcode at the right end of the read. Since we normally expect the barcode on the left, we use that to assign the read (and still trim both barcodes).

**In short:** if you want to be conservative, do not use `--maximize`. If a small number of potential misassignments is acceptable, `--maximize` can help retain more reads.

## Custom experiment

> [!IMPORTANT]
> If trimming fails/crashes you likely have many samples and need to increase ulimit `ulimit -n 65535`

We first create a Fasta, or multiple Fastas containing your queries depending on whether 
you have a single-end or dual-end experiment.


### Creating a query Fasta
If you have your own barcodes/primers/adapters, create FASTA files containing the full expected sequences. 
Note that each FASTA contains all sequences that <u>share the same prefix/suffix</u>, and are **unique**. 

The Fasta format should be as follows:
```text
>NB01
<left_flank_sequence><bar1_sequence><right_flank_sequence>
>NB02
<left_flank_sequence><bar2_sequence><right_flank_sequence>
...
```

Example (adapter + barcode + primer):

```
>NB01
TCGTTCAGTTACGTATTGCTCACAAAGACACCGACAACTTTCTTAGRGTTYGATYATGGCTCAG
>NB02
TCGTTCAGTTACGTATTGCTACAGACGACTACAAACGGAATCGAAGRGTTYGATYATGGCTCAG
```

Here:
```
TCGTTCAGTTACGTATTGCT CACAAAGACACCGACAACTTTCTT AGRGTTYGATYATGGCTCAG
    adapter               barcode                 primer
```

Sarracenia extracts the shared prefix (adapter) and shared suffix (primer) as flanks — only the barcode region should differ between FASTA entries.


### Single end
See above to create a fasta file, say `left.fasta` for single-end then we can run `annotate`:

```bash
sarracenia annotate -q left.fasta -b Ftag -i reads.fastq -o anno.tsv -t 10
```
After that you can follow the same `inspect → filter → trim` [steps described above](#inspect).


### Dual end
For dual-end, we create two FASTAs, `left.fasta` and `right.fasta` for example, then we run:

```bash
sarracenia annotate -q left.fasta,right.fasta -b Ftag,Rtag -i reads.fastq -o anno.tsv -t 10
```
Note: there must be **no spaces** between the comma-separated file list and the tag list: `-q left.fasta,right.fasta -b Ftag,Rtag`.


In case you have concatenated reads in your `inspect` also see [concat reads](#how-to-handle-concat-reads)





## Custom experiment with mixed sequences
Often we combine multiple samples together which share the same flanks, for example all rapid barcoding. 
But we could also combine completely different experiments, with different primers for example.
How do we demux that?

First go through "[custom experiment](#custom-experiment)" setup, then continue reading here how we adjust the query files.

That's quite simple. We create fasta files for all possible groups. 
Lets say we have:

**group1**
- group1_left.fasta: adapter1-barcode-primer1
- group1_right.fasta: primer1-barcode-adapter1

**group2**:
- group2_left.fasta: adapter2-barcode-primer2
- group2_right.fasta: primer2-barcode-adapter2

Then we <u>make sure our labels in the fasta have some unique substring</u>, like `group1` and `group2`:
for example:

`group1_left.fasta`:
```
>group1_bar1
AACGACA...
```
`group2_left.fasta`:
```
>group2_bar1
AGGGCAC...
```

We then run annotate like:

```
sarracenia annotate -q group1_left.fasta,group1_right.fasta,group2_left.fasta,group2_right.fasta -b Ftag,Rtag,Ftag,Rtag -i reads.fastq -o anno.tsv -t 10
```

**Note**: While you could run annotate separately for each query file, it’s generally better to combine all into a single annotate run (as mentioned here). This way, they are competing with one another.


Then in the filtering step we can create filtered files for each group:
(just check `inspect` first to see your patterns)


`group1_filters.txt`:
```
Ftag[fw, ~group1, @left(0..250), >>]__Rtag[<<, rc, ~group1, @right(0..250)]
```

`group2_filters.txt`:
```
Ftag[fw, ~group2, @left(0..250), >>]__Rtag[<<, rc, ~group2, @right(0..250)]
```

And pull them out!
```
sarracenia filter -i anno.tsv -f group1_filters.txt -o group1_reads.tsv
sarracenia filter -i anno.tsv -f group2_filters.txt -o group2_reads.tsv
```
and then we can just trim them to separate files:

```
sarracenia trim -i group1_reads.tsv -r reads.fastq -o group1_trimmed
sarracenia trim -i group2_reads.tsv -r reads.fastq -o group2_trimmed
```





## Output columns (annotate & filter)

- `read_id`: read identifier as in the provided FASTQ
- `read_len`: length of read in bp
- `rel_dist_to_end`: relative distance to the read end. `>= 0` means X bases from the **left** end; a negative value means X bases from the **right** end (e.g. `-10` is 10 bp from the right).
- `read_start_bar`: start position of the barcode match in the read
- `read_end_bar`: end position of the barcode match in the read
- `read_start_flank`: start position of the flank match in the read
- `read_end_flank`: end position of the flank match in the read
- `bar_start`, `bar_end`: coordinates where the barcode was aligned (previous partial-barcode matches are disabled; full barcode length is expected)
- `match_type`:
  - `Ftag`: forward barcode + flank matched
  - `Fflank`: forward flank matched (barcode undetectable)
  - `Rtag`: rear barcode + flank matched
  - `Rflank`: rear flank matched (barcode undetectable)

  *Note*: for some kits (e.g., rapid barcoding) there's effectively a single barcode; we still call this `Ftag`. For native dual barcoding both ends may use the same barcode set; orientation (forward/reverse complement) is available in the `strand` column.

- `flank_cost`: number of edits in the flank sequence (excluding the barcode)
- `barcode_cost`: number of edits in the barcode region
- `label`: label from your FASTA (or preset kit) — e.g., `BC14`, `RBK60`, etc.
- `strand`: orientation of the match (`fw` or `rc`)
- `cuts`: empty after `annotate`; populated after `filter` to inform `trim` where to cut





## Patterns

Patterns describe how elements are combined in a read. Single elements have the form:

```
<tag>[<orientation>, <label>, <relative position>, <optional cut specifier>]
```

Multiple elements are combined with `__` (double underscore):

```
<tag>[...]__<tag>[...]
```

Fields:
- `tag`: `Ftag`, `Fflank`, `Rtag`, `Rflank`, whether a sequence is `Ftag` or `Rtag` is user specified (or within the kit), see [custom experiment](#custom-experiment), the `Fflank,Rflank` are the "incomplete" forms where the barcode was undetectable.
- `orientation`: `fw` or `rc`.
- `label`: exact label (e.g. `NB01`), `*` for any label, or `~substring` to match headers containing `substring` (for an example see [custom exp. mixing](#custom-experiment-with-mixed-sequences)).
- `relative position`: e.g. `@left(0..250)`: to left side of read, `@right(0..250)`: to right side of read, `@prev_left(0..250)`: relative to the <u>previous</u> element.
- `cut specifier` (optional): `>>` (cut after this element) or `<<` (cut before this element).

Examples:

```
Ftag[fw, *, @left(0..250), >>]
Ftag[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]
Ftag[fw, *, @left(0..100), >>]__Rtag[<<, fw, *, @prev_left(1500..1700)]
```
(for examples of concat reads see [concat reads](#how-to-handle-concat-reads) section)

These let you express typical cases such as single-barcode-left, left-and-right barcodes, or expected amplicon sizes via `@prev_left`.




### How to handle concat reads 
It is possible that you have concat reads like:

```
Ftag[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @pev_left(1500..1750)]__Ftag[fw, *, @prev_left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]
```
We always read the pattern from <u> left to right </u>, then we can deduce the read matches looked like this:
```
[Ftag,fw][Ftag,rc][Ftag, fw][Ftag, rc]
```
If we want to cut more than once we have to use **cut group identifiers**, we do this by placing a number <u> after </u> the `<<` or `>>`:

```
Ftag[fw, *, @left(0..250), >>1]__Ftag[<<1, rc, *, @pev_left(1500..1750)]__Ftag[fw, *, @prev_left(0..250), >>2]__Ftag[<<2, rc, *, @right(0..250)]
```
Note the `>>1, <<1, >>2, <<2`, now Sarracenia knows exactly where you want the reads to be cut.

If the labels were: `NB01NB01-NB02-NB02`, Sarracenia will write the first read to the `NB01` folder and the second to `NB02`. 

#### Getting *all* concats

If concats are a big part of your read set (ideally they are not) their patterns will show up high in `inspect` anyway, and copying just those will 
probably be enough. If you really want *all* the concat reads out you could use the `-o` in `inspect` to write the patterns
per read. Then, get all unique patterns from there and use a regex or a python script to insert the `>>1` and `<<1`, you can just dump all of them 
to `filters.txt` and use that in `sarracenia filter`. 




## Paper evals
Since this involves substantial amount of extra code we moved these to the [paper-evals](https://github.com/katalyxbio/sarracenia-evals) repo. All information on how to reproduce results and set up environments to do so can be found there. 



## Notes & tips

- Terminal output is plain text with no colour or ANSI art. Informational lines go to
  stdout; each stage's progress line goes to stderr, so `sarracenia ... > run.log` keeps
  the two apart. On a terminal the progress line rewrites itself in place; when stderr is
  redirected it is appended every 10 seconds instead, so logs stay short. Pass the global
  `--quiet` flag to print nothing but errors.
- Keep an eye on `Fflank` matches — these are often lower-confidence and may indicate reads with poor barcode sequence quality.
- Start with conservative filters to see how many reads match expected patterns, then relax thresholds if necessary. 
- When experimenting with custom FASTAs, keep queries clean (shared flanks, differing barcode region only).

## License

Sarracenia is licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0),
see [LICENSE](LICENSE). You may use, modify and share it for **noncommercial purposes** —
which the license defines as personal use, and use by charitable organizations, educational
institutions, public research organizations, public safety or health organizations,
environmental protection organizations, and government institutions.

Using Sarracenia for a commercial purpose requires a separate license from the copyright
holder. Note that this is a source-available, not an open-source, license: it is not
OSI-approved.
