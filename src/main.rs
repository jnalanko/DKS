#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter}, path::PathBuf, sync::Arc};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
use indicatif::HumanBytes;
use io::LazyFileSeqStream;
use jseqio::reader::DynamicFastXReader;
use sbwt::{BitPackedKmerSorting, BitPackedKmerSortingMem, LcsArray, SbwtConstructionAlgorithm, SbwtIndex, SbwtIndexVariant, SeqStream, SubsetMatrix};
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

        #[arg(help = "Directory for temporary files. This enables disk-based construction.", short = 'd', long = "temp-dir")]
        temp_dir: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", short = 'f', long = "forward-only")]
        forward_only: bool,

        #[arg(short, required = true)]
        k: usize,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },

    #[command(arg_required_else_help = true)]
    Lookup {
        #[arg(help = "A file with one fasta/fastq filename per line", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,

        //#[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        //n_threads: usize,
    },
}


fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Cli::parse();
    match args.command {
        Subcommands::Build { input: input_fof, output: out_path, temp_dir, k, n_threads, forward_only} => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof).unwrap()).lines().map(|f| PathBuf::from(f.unwrap())).collect();
            let mut out = BufWriter::new(File::create(out_path.clone()).unwrap());

            let all_input_seqs = io::ChainedInputStream::new(input_paths.clone());
            let add_rev_comps = !forward_only;

            let (sbwt, lcs) = if let Some(td) = temp_dir {
                // Use disk-based construction
                sbwt::SbwtIndexBuilder::new()
                    .add_rev_comp(add_rev_comps)
                    .k(k)
                    .build_lcs(true)
                    .n_threads(n_threads)
                    .precalc_length(8)
                    .algorithm(BitPackedKmerSorting::new().dedup_batches(false).temp_dir(&td))
                .run(all_input_seqs)
            } else {
                // Use in-memory construction
                sbwt::SbwtIndexBuilder::new()
                    .add_rev_comp(add_rev_comps)
                    .k(k)
                    .build_lcs(true)
                    .n_threads(n_threads)
                    .precalc_length(8)
                    .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
                .run(all_input_seqs)
            };
            let lcs = lcs.unwrap(); // Ok because of build_lcs(true)

            let individual_streams = input_paths.iter().map(|p| LazyFileSeqStream::new(p.clone(), add_rev_comps)).collect();
            log::info!("Marking colors");
            let index = SingleColoredKmers::new(sbwt, lcs, individual_streams, n_threads);

            log::info!("Writing to {}", out_path.display());
            index.serialize(&mut out);
            let out_size = std::fs::metadata(out_path).unwrap().len() as f64;
            log::info!("Index size on disk: {}" , human_bytes::human_bytes(out_size));
        },

        Subcommands::Lookup{query: query_path, index: index_path} => {
            eprintln!("Loading the index");
            let mut index_input = BufReader::new(File::open(index_path).unwrap());

            let index_loading_start = std::time::Instant::now();
            let index = SingleColoredKmers::load(&mut index_input);
            eprintln!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());
            let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
            let mut color_hit_counts = vec![0_usize; index.n_colors()];
            let mut total_kmers_queried = 0_usize;
            let query_start = std::time::Instant::now();
            while let Some(rec) = reader.read_next().unwrap() {
                for color in index.lookup_kmers(rec.seq) {
                    if let Some(color) = color {
                        color_hit_counts[color] += 1;
                    }
                    total_kmers_queried += 1;
                }
            }

            let query_duration = query_start.elapsed();
            for color in 0..index.n_colors() {
                let hits = color_hit_counts[color];
                println!("Color {}: {} hits ({:.2}%)", color, hits, hits as f64 / total_kmers_queried as f64 * 100.0);
            }
            
            eprintln!("{} k-mers queried in {} seconds (excluding index loading time)", total_kmers_queried, query_duration.as_secs());
            eprintln!("{:.2}% of query k-mers found", color_hit_counts.iter().sum::<usize>() as f64 / total_kmers_queried as f64 * 100.0);
            eprintln!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / total_kmers_queried as f64);
        }
    } 
}