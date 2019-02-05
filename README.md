# bamcount

BigWig and BAM related utilities.

## Build dependencies

* [htslib](http://www.htslib.org)
    * See `get_htslib.sh` for a script that gets a recent version and compiles it with minimal dependencies
* [libBigWig](https://github.com/dpryan79/libBigWig)
    * See `get_libBigWig.sh` for a script that gets a recent version anf compile it

## Compiling

From root directory:

```
mkdir -p build && cd build && cmake .. && make
```

## Usage

```
bamcount --coverage /path/to/bamfile --threads <num_threads> --no-head --coverage --bigwig <sample_name> --auc --min-unique-qual <min_qual> --annotation <annotated_intervals.bed> <sample_name> --frag-dist <sample_name> --alts <sample_name>
```

Concrete example command for sample `SRR1258218` (NA12878 Illumina RNA-seq):

```
bamcount ./SRR1258218.sorted.bam --threads 4 --no-head --coverage --bigwig SRR1258218 --auc SRR1258218 --min-unique-qual 10 --annotation ./exons.bed SRR1258218 --frag-dist SRR1258218 --alts SRR1258218
```

## Subcommands

For any and all subcommands below, if run together, `bamcount` will do only one pass through the BAM file.
While any given subcommand may not be particularly fast on its own, doing them all together can save time.

Subcommand `--coverage` is the only subcommand that will output a BigWig file (currently).
However, if along with `--coverage`, `--min-unique-qual` and `--bigwig` are specified the "unique" coverage will also be written as a BigWig file.

### `bamcount --auc <output_file_prefix>`

Reports area-under-coverage across all bases (one large sum of overlapping reads, per-base).
This will also report additional counts for:
 * `min-unique-qual` only for reads with MAPQ >= to this setting
 * `--annotation`: only for bases in the annotated regions
 
This computes the coverage (same as `--coverage`) under the hood, but won't output it unless `--coverage` is also passed in.

### `bamcount --coverage`

Generates per-base counts of overlapping reads across all of the genome.  
Typically this is used to produce a BigWig.

### `bamcount --coverage --min-unique-qual <qual_value>`

In addition to producing coverage output for all reads, will also produce coverage output only for reads which have a mapping quality (MAPQ) >= <qual_value> (typically set to `10`.  

### `bamcount --coverage --annotation <annotated_file.bed> <output_file_prefix>`

In addition to reporting per-base coverage, this will also sum the per-base coverage within annotated regions submitted as a BED file.
If `--min-unique-qual` is submitted, this will produce a second set of sums for the "unique" reads that pass this filter.
 
### `bamcount --coverage --double-count`

By default, `bamcount --coverage` will not double count coverage where paired-end reads overlap (same as `mosdepth`'s default).
However, double counting can be allowed with this option, which may result in faster running times.

### `bamcount --coverage --bigwig <output_file_prefix>`

Outputs coverage vectors as BigWig file(s) (including for `--min-unique-qual` option).

### `bamcount --frag-dist <output_file_prefix>`

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

### `bamcount --alts <output_file_prefix>`

Outputs information about non-reference-matching portions of reads.
Output is comma separated with 4 fields:

| Pos   |                                            Descrtiption |
|-------|---------------------------------------------------------|
| 1     | Reference record ID                                     |
| 2     | POS field (1-based offset of leftmost aligned ref base) |
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

