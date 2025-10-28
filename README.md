# Disjoint k-mer sets

DKS (Disjoint Kmer Sets) is a tool for indexing sets of k-mers that do not intersect. Each input k-mer set is assigned a unique *color* that indicates which set the k-mer is from. The query takes an index and a set of query sequences, and reports the color of each k-mer in the query.

The indexing is based on the [Spectral Burrows-Wheeler Transform](https://docs.rs/sbwt/latest/sbwt).

## Installation

Precompiled x86-64 Linux binaries are available at [releases](https://github.com/jnalanko/DKS/releases).

To build from source, follow these steps:

* Step 1: [Install the latest Rust compiler](https://rust-lang.org/tools/install), if not already installed. This does not require root privileges.
* Step 2: Clone and enter the repository: `git clone https://github.com/jnalanko/DKS --recursive; cd DKS`.
* Step 3: Compile and install by running `cargo install --path .` in this directory. 

Now you should be able to run DKS by typing `dks` on the command line.

## Quick start

**Indexing**. This repository contains a small example dataset to demonstrate the usage of the tool. In the subdirectory `example`, there are files `A.fasta`, `B.fasta` and `C.fasta`. The input for indexing is a file with one filename for each line, like so:

```
./example/A.fasta
./example/B.fasta
./example/C.fasta
```

This file is at `example/file_of_files.txt`. The files are assigned to colors 0,1,2 in the same order as they appear in the list. All k-mers in a file are colored with the color of the file. If the same k-mer appears in multiple files, it will be given a special color indicating that it has multiple colors, but the index does not store which colors those are.

To index this data with k = 5, using 8 threads and the directory `./temp` for temporary working space, run: 

```bash
dks build -k 5 -i example/file_of_files.txt -o example/index.dks --external-memory ./temp -t 8
```

This will save the index to `example/index.dks`. Reverse complements of k-mers are automatically indexed and receive the same color as the forward k-mer. In a real use case, you might want to use something like k = 31. **Too slow?** Remove the option `--external-memory ./temp` to run entirely in RAM. See [benchmarks](#performance) for how this affects time and space on a human genome. See also [best practices](#best-practices) for more advice on how to maximize performance.

**Lookup**. To query all k-mers in the file `example/query.fasta` against the index built above, using 8 threads, writing the output to `example/out.tsv`, run:

```bash
dks lookup -q example/query.fasta -i example/index.dks -t 8 > example/out.tsv
```

This will write the following to `example/out.tsv`:

```
seq_rank	from_kmer	to_kmer	color
0	0	0	*
0	1	2	0
0	3	3	*
0	4	4	1
0	7	10	2
0	11	11	*
0	12	14	2
1	0	0	*
1	1	7	1
1	12	13	0
1	14	14	*
1	19	21	2
```

The four columns are:

* `seq_rank`: Zero-based sequence rank in the input file
* `from_kmer`: Zero-based starting position of the first k-mer of a k-mer range
* `to_kmer`: Zero-based starting position of the last k-mer of a k-mer range (inclusive)
* `color`: Color of the k-mers in the range, or `*` if the k-mer has multiple colors

So in our example, for the first query sequence (sequence 0), k-mers 0-0 have multiple colors, k-mers 1-2 have color 0, k-mers 3-3 have multiple colors, and so on.

## Performance

### Indexing and query results on the human genome

In this benchmark, all chromosomes of the human reference genome GRCh38.p14 were indexed, such that each chromosome received a distinct color (24 distinct colors total). Then, a full human genome assembly was queried against the index. 

![Benchmark plots](benchmark/benchmarks_combined.png)

Theoretical ideal parallel speedups are shown with dashed lines. The size of the index on disk was 8.8 GiB for k = 31 and 10.3 GiB for k = 63. The peak memory during query was 9.0 GiB for k = 31 and 10.6 GiB for k = 63. The experiment was run on a SLURM cluster on a node with an AMD EPYC 7452 32-Core processor clocked at 2350 MHz with cache sizes 32K, 512K and 16384K, and 500GB of RAM (detailed specs unknown).

### Indexing memory estimate

When run without `--external-memory`, the memory is roughly 16n ⌈k/32⌉ bytes, where n is the number of nucleotides in the input. For example, when indexing the 31-mers in all chromosomes of the human genome (3.3 billion nucleotides), the predicted space is 49.2 GiB, and the actual space is 51.2 GiB. When run with external memory enabled, this amount of data is written and read from the disk instead, and the peak memory is only 13.4 GiB.

In pathological cases such as when the k-mers do not overlap at all, the space can be larger. 

## Best practices 

* When running with `--external-memory`, make sure that the working directory given as the parameter is set to a location with fast sustained writes and reads, preferably a fast SSD drive.
* An underlying assumption in the tool is that the input k-mers come from longer contiguous sequences. This is exploited to improve running time and space. For best performance, give the inputs as long sequences of overlapping k-mers instead of giving every k-mer separately.

## Detailed instructions 

### Build

```
Usage: dks build [OPTIONS] -k <K> --input <INPUT> --output <OUTPUT>

Options:
  -k <K>                            
  -i, --input <INPUT>               A file with one fasta/fastq filename per line
  -o, --output <OUTPUT>             Output filename
      --external-memory <TEMP_DIR>  Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.
  -f, --forward-only                Do not add reverse complemented k-mers
  -t, --n-threads <N_THREADS>       Number of parallel threads [default: 4]
  -s, --sbwt-path <SBWT_PATH>       Optional: a precomputed SBWT file of the input k-mers.
  -l, --lcs-path <LCS_PATH>         Optional: a precomputed LCS file of the optional SBWT file.
  -h, --help                        Print help
```

### Lookup

```
Usage: dks lookup [OPTIONS] --query <QUERY> --index <INDEX>

Options:
  -q, --query <QUERY>          A file with one fasta/fastq filename per line
  -i, --index <INDEX>          Path to the index file
  -t, --n-threads <N_THREADS>  Number of parallel threads [default: 4]
  -h, --help                   Print help
```

## For developers: Building portable Linux binaries for release

For Linux, install the build toolchain for target x86_64-unknown-linux-musl:

```bash
rustup target add x86_64-unknown-linux-musl
```

Compile with:

```bash
cargo build --release --target x86_64-unknown-linux-musl
```

