#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter}, path::PathBuf, sync::Arc};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
use io::LazyFileSeqStream;
use sbwt::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SbwtIndexVariant, SubsetMatrix};
use single_colored_kmers::SingleColoredKmers;

mod single_colored_kmers;
mod io;


#[derive(Parser)]
#[command(arg_required_else_help = true)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Subcommands,
}

#[derive(Subcommand)]
pub enum Subcommands {
    #[command(arg_required_else_help = true)]
    Build {
        #[arg(help = "A file with one fasta/fastq filename per line", short, long, required = true)]
        input: PathBuf,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Directory for temporary files", short = 'd', long = "temp-dir", required = true)]
        temp_dir: PathBuf,

        #[arg(short, required = true)]
        k: usize,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },
}


fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads} => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();

            let all_input_seqs = io::ChainedInputStream::new(input_paths.clone());
            let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
                .add_rev_comp(true)
                .k(k)
                .build_lcs(true)
                .n_threads(n_threads)
                .precalc_length(8)
                .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
            .run(all_input_seqs);
            let lcs = lcs.unwrap(); // Ok because of build_lcs(true)

            let individual_streams = input_paths.iter().map(|p| LazyFileSeqStream::new(p.clone())).collect();
            let index = SingleColoredKmers::new(sbwt, lcs, individual_streams);
            // TODO: add reverse complements
            //index.serialize();
        },
    } 
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_to_classes(){
        let mut vec = vec![
            bitvec![1, 1, 1],
            bitvec![1, 0, 1],
            bitvec![0, 1, 0],
            bitvec![1, 0, 1],
            bitvec![1, 0, 0],
            bitvec![0, 1, 0],
            bitvec![0, 1, 0],
            bitvec![1, 0, 1]];
        let counts = reduce_to_classes(&mut vec);
        assert_eq!(vec, vec![
            bitvec![0, 1, 0],
            bitvec![1, 0, 0],
            bitvec![1, 0, 1],
            bitvec![1, 1, 1]]);
        assert_eq!(counts, vec![3, 1, 3, 1]);
    }
}
