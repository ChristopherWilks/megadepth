# megadepth

![build](https://github.com/ChristopherWilks/megadepth/workflows/build/badge.svg)

BigWig and BAM related utilities.

We strongly recommend use of one of the pre-compiled binaries for x86_64 linux systems:

* [dynamically linked binary with HTSlib, libBigWig, libcurl, libdeflate, & zlib statically linked](https://github.com/ChristopherWilks/megadepth/releases/download/1.0.5/megadepth)

NOTE: the statically linked binary does not support remote BigWig/BAM processing due to the difficulties in linking a static libcurl, but may still be useful for those who want to do local processing on systems where the dynamic binary doesn't work.

There is also a Docker image that can be used to run `megadepth`:

https://quay.io/repository/broadsword/megadepth?tab=tags

You'll probably want to map in a directory on the host system into the container via the `-v` option so you can pass an annotation file in and get output back:

```
docker run -v `pwd`:/data <image_id> </data/path_or_URL/to/input/BAM_or_BigWig> --annotation /data/<annotation>.bed --prefix /data/output_file_prefix
```

Currently, `libcurl` throws a warning about version information, this can be ignored.

Finally, if none of those options work, the build instructions are at the end of this README.

[Releases prior to 1.0.2 used the previous name "bamcount"]

## Usage

For any remote file processing, either BAM or BigWigs, you *must* use the `--prefix <output_file_prefix>` option.

### BigWig Processing
```
megadepth /path/to/bigwigfile --annotation <annotated_intervals.bed> --op <operation_over_annotated_intervals>
```

Concrete example command for sample `SRR1258218` (NA12878 Illumina RNA-seq), this will produce 1) means for the intervals listed in `exons.bed` and 2) the total annotated AUC (output `STDOUT`):
```
megadepth SRR1258218.bw --annotation exons.bed --op mean --auc
```

Or if you only want the AUC for the whole BigWig:
```
megadepth SRR1258218.bw
```

### BAM processing
While megadepth doesn't require a BAM index file (typically `<prefix>.bam.idx`) to run, it *does* require that the input BAM be sorted by chromosome at least.  This is because megadepth allocates a per-base counts array across the entirety of the current chromosome before processing the alignments from that chromosome.  If reads alignments are not grouped by chromosome in the BAM, undefined behavior will occur including massive slow downs and/or memory allocations.

```
megadepth /path/to/bamfile --threads <num_threads> --bigwig --auc --annotation <annotated_intervals.bed> --prefix <output_file_prefix>
```

Concrete example command for sample `SRR1258218` (NA12878 Illumina RNA-seq):

```
megadepth SRR1258218.sorted.bam --threads 4 --bigwig --auc --annotation exons.bed --prefix SRR1258218
```

## BAM Processing Subcommands

For any and all subcommands below, if run together, `megadepth` will do only one pass through the BAM file.
While any given subcommand may not be particularly fast on its own, doing them all together can save time.

Subcommand `--bigwig` is the only subcommand that will output a BigWig file with the suffix `.all.bw`.
If `--min-unique-qual` and `--bigwig` are specified the "unique" coverage will also be written to a separate BigWig file with the suffix `.unique.bw`.

Also, `--bigwig` will not work on Windows, megadepth as of release 1.0.5 will simply skip writing a BigWig if this option is passed in with the Windows build, but will process other options which still make sense (e.g. `--auc`).

### `megadepth /path/to/bamfile --auc`

Reports area-under-coverage across all bases (one large sum of overlapping reads, per-base).
This will also report additional counts for:
 * `min-unique-qual` only for reads with MAPQ >= to this setting
 * `--annotation`: only for bases in the annotated regions
 
This computes the coverage (same as `--coverage`) under the hood, but won't output it unless `--coverage` is also passed in.

Will default to reporting to `STDOUT` unless `--no-auc-stdout` is passed in.

### `megadepth /path/to/bamfile --coverage`

Generates per-base counts of overlapping reads across all of the genome.  
Typically this is used to produce a BigWig, but can be used w/o the `--bigwig` option to just output TSVs

Will default to reporting to `STDOUT` unless `--no-coverage-stdout` is passed in.

### `megadepth /path/to/bamfile --coverage --annotation <annotated_file.bed> --no-coverage-stdout --no-annotation-stdout`

In addition to reporting per-base coverage, this will also sum the per-base coverage within annotated regions submitted as a BED file.

The annotation BED file does not need to be sorted in any particular way.

megadepth will output the summed coverages for the annotation in contiguous blocks per chromosome.

This will be the same order as the BED file *if* coordinates from the same chromosome are contiguous in the BED file (typically they are).

It's best to use the command as given, otherwise both the coverage and annotated coverage output will be reported intermingled to `STDOUT`.

Will default to reporting to `STDOUT` unless `--no-annotation-stdout` is passed in.

Also, this no longer automatically reports the AUC, you'll also need to pass in `--auc` if you want that as well.
 
### `megadepth /path/to/bamfile --coverage --double-count`

By default, `megadepth --coverage` will not double count coverage where paired-end reads overlap (same as `mosdepth`'s default).
However, double counting can be allowed with this option, which may result in faster running times.

### `megadepth /path/to/bamfile --bigwig`

Outputs coverage vectors as BigWig file(s) (including for `--min-unique-qual` option).

### `megadepth /path/to/bamfile --frag-dist`

Outputs fragment length distribution adjusting for intron lengths.

Mean, mode statistics are reported at the end of the output with string tag `STATS`.

This uses the absolute value of the `TLEN` field but uses additional filters similar to [csaw](https://github.com/LTLA/csaw)'s fragment length calculation.

The following alignments are filtered out:

 * secondary
 * supplementary
 * not paired
 * unmapped
 * mate unmapped
 * discordant (mates not on same chromosome/reference)

Further, read mates must be on forward/reverse strands and the forward mate must not be downstream of the reverse mate.

Intron length(s) in the paired alignments are also subtracted from the `TLEN` field except where the `TLEN` field is smaller than the combined length of the introns, in which case the `TLEN` is reported as is.

These numbers should be taken as an estimation of the fragment length distribtion.

Reports to a file with suffix `.frags.tsv`.

### `megadepth /path/to/bamfile --alts`

Outputs information about non-reference-matching portions of reads.
Output is comma separated with 4 fields:

| Pos   |                                            Descrtiption |
|-------|---------------------------------------------------------|
| 1     | Reference/chromosome ID in the BAM file (integer)       |
| 2     | POS field (0-based offset of leftmost aligned ref base) |
| 3     | Operation label (see table below)                       |
| 4     | Extra info (see table below)                            |


These could be of a few types, summarized in this table.  All of
these are available when the `MD:Z` extra flag is present.  If not
present, only the ones with "Yes" in the "No `MD:Z`" column are
reported.

| Label | Type       |          Extra info  | No `MD:Z`  |
|-------|------------|----------------------|------------|
| `X`   | Mismatch   |            Read base |         No |
| `D`   | Deletion   |      # deleted bases |        Yes |
| `I`   | Insertion  |  Inserted read bases |        Yes |
| `S`   | Soft clip  |   Soft ckipped bases |        Yes |
| `H`   | Hard clip  |            (nothing) |        Yes |
| `P`   | Padding    |            (nothing) |        Yes |

See the usage message for options, which can selectively disable some
of the outputs listed above.  E.g. the soft-clipping outputs can be
very large, so they're not printed unless `--include-softclip` is
specified.

Reports to a file with suffix `.alts.tsv`.

### `megadepth /path/to/bamfile --alts --include-softclip`

In addition to the alternate base output, this reports the bases
that were softclipped at the ends (start/end) of the read.
These are bases which are left in the sequence but don't align.

The softclipped bases themselves are printed to the file named with the
prefix passed into the `--alts` option.  The total number of sofclipped
bases and the total number of bases from the query sequences of alignments
that that aren't unmapped or secondary are reported to the file named
with the prefix passed to `--include-softclip`.

Warning: using this option w/o modifiers (e.g. `--only-polya`) 
could blow up the `--alts` output size as the full softclipped
sequence is printed in the 4th column in the table above ("Extra info").

Reports to a file with suffix `.softclip.tsv` in addition to the `--alts` file.

### `megadepth /path/to/bamfile --alts --include-softclip --only-polya`

If reporting softclipped bases, this option will limit the report to only
those bases that have the following:

* Count of bases in the sofclip (column 4 below) has to be >= 3
* % of base (A/T) of softclipped bases for an alignment >= 80%

No other sofclipped bases are reported.

Output is comma separated with 7 fields:

| Pos   |                                                                      Description|
|-------|----------------------------------------------------------------------------------|
| 1     | Reference/chromosome ID in the BAM file (integer)                                |
| 2     | POS field (0-based ref offset of either leftmost or rightmost aligned base)      |
| 3     | Operation label (always "S")                                                     |
| 4     | Number of bases in the softclip (run length)                                     |
| 5     | Direction to move from POS ('+' for end of alignment, '-' for start of alignment)|
| 6     | Base (A/T)                                                                       |
| 7     | Count of the base in column 6                                                    |

### `megadepth /path/to/bamfile --junctions`

Extract locally co-occurring junctions from BAM.

This does not extract all potential junctions, only those for which a read (or read pair) had >= 2 junctions.

In a paired context, there must be at least 2 junctions across the 2 read mates to be output.

Output is tab separated with 6-12 fields (the last 6 fields are for a 2nd mate if applicable):

| Pos    |                                                                            Description|
|--------|---------------------------------------------------------------------------------------|
| 1      | Reference/chromosome ID in the BAM file (integer)                                     |
| 2      | POS field (1-based ref offset of either leftmost base)                                |
| 3      | Mapping strand (0 forward, 1 reverse)                                                 |
| 4      | Insert length (0 if not paired)                                                       |
| 5      | Cigar string (useful for determining anchor lengths)                                  |
| 6      | List of junction coordinates (comma-delimited)                                        |
| 7*     | Mate reference record ID                                                              |
| 8*     | Mate POS field (1-based ref offset of either leftmost base)                           |
| 9*     | Mate mapping strand (0 forward, 1 reverse)                                            |
| 10*    | Mate insert length (0 if not paired)                                                  |
| 11*    | Mate cigar string (useful for determining anchor lengths)                             |
| 12*    | Mate list of junction coordinates (comma-delimited)                                   |

\*optional, output if a 2nd mate is present and has the required number of junctions.

If you get a core dump when running on longer reads (e.g. BAM's produced by PacBio/Oxford Nanopore sequencing),
then try adding the argument `--long-reads` as it will enlarge the buffer used to store the output junction string.

This enables megadepth to have a better chance of handling really long CIGAR strings.

Reports to a file with suffix `.jxs.tsv`.

## Build dependencies

* [htslib](http://www.htslib.org)
    * See `get_htslib.sh` for a script that gets a recent version and compiles it with minimal dependencies
* [libBigWig](https://github.com/dpryan79/libBigWig)
    * See `get_libBigWig.sh` for a script that gets a recent version and compiles it
* zlib static library [only if building a static binary]
    * See `get_zlib.sh` for a script that gets a recent version and compiles the static library

## Building

Run `build_no_container.sh` with one of three options:

* `megadepth_dynamic` (default)

Builds a fully dynamic binary, requires that libraries for `htslib` & `libBigWig` be available in the target environment

* `megadepth_statlib`

Builds a partially dynamic binary, but with `htslib` and `libBigWig` statically linked, still requires that libcurl and zlib be present in the target environment

* `megadepth_static`

Builds a fully static binary, w/o remote BigWig processing support (due to no libcurl)


