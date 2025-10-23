#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{cmp::{max, Reverse}, fs::File, io::{BufRead, BufReader, BufWriter}, ops::Range, path::{Path, PathBuf}, sync::Arc, thread::ScopedJoinHandle};
use bitvec::prelude::*;
use clap::{builder::PossibleValuesParser, Parser, Subcommand};
use indicatif::HumanBytes;
use io::LazyFileSeqStream;
use jseqio::{reader::DynamicFastXReader, seq_db::SeqDB};
use sbwt::{BitPackedKmerSorting, BitPackedKmerSortingMem, ContractLeft, ExtendRight, LcsArray, SbwtConstructionAlgorithm, SbwtIndex, SbwtIndexVariant, SeqStream, StreamingIndex, SubsetMatrix};
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

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
    },
}

struct QueryBatch {
    seqs: SeqDB,

    // The sequences in `seqs` are subsequence of some longer sequences S_1, ... S_n.
    // The first k-mer in the first sequence in `seqs` starts at `start_kmer_offset` in S_{seq_id_range.start}
    seq_id_range: Range<usize>,
    start_kmer_offset: usize,
    end_kmer_offset: usize, // One past the starting point of the last k-mer

    result: Vec<Option<usize>>, // Color ids of the query k-mers
}

impl QueryBatch {
    fn run(&mut self, index: &SingleColoredKmers) {
        let k = index.k();
        let total_query_kmers = self.seqs.iter().fold(0_usize, |acc, rec| 
            acc + kmers_in_n(k, rec.seq.len()) 
        );
        let mut color_ids = Vec::<Option::<usize>>::with_capacity(total_query_kmers as usize);

        for rec in self.seqs.iter() {
            for color in index.lookup_kmers(rec.seq) {
                if let Some(color) = color {
                    color_ids.push(Some(color));
                } else {
                    color_ids.push(None);
                }
            }
        }

        assert_eq!(color_ids.len(), total_query_kmers as usize);
        self.result = color_ids;

    }
}

impl PartialEq for QueryBatch {
    fn eq(&self, other: &Self) -> bool {
        self.seq_id_range == other.seq_id_range && self.start_kmer_offset == other.start_kmer_offset
    }
}
impl Eq for QueryBatch {}

impl PartialOrd for QueryBatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other)) // Using the total order
    }
}

impl Ord for QueryBatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Primary key: seq id range start, secondary key: kmer start
        (self.seq_id_range.start, self.start_kmer_offset).cmp(&(other.seq_id_range.start, other.start_kmer_offset))
    }
}

fn kmers_in_n(k: usize, n: usize) -> usize {
    max(0, n as isize - k as isize + 1) as usize
}

fn output_thread(query_results: crossbeam::channel::Receiver<QueryBatch>) {
    let cur_seq_id = 0_usize;
    let cur_kmer_offset = 0_usize; 

    let mut batch_buffer = std::collections::BinaryHeap::<Reverse<QueryBatch>>::new(); // Reverse makes it a min heap

    while let Ok(batch) = query_results.recv() {
        batch_buffer.push(Reverse(batch)); // Reverse makes this a min heap

        loop {
            let min_batch = batch_buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.seq_id_range.start == cur_seq_id && min_batch.start_kmer_offset == cur_kmer_offset {
                    for x in min_batch.result.iter() {
                        println!("{:?}", x); // TODO: better printing
                    }
                    batch_buffer.pop();
                } else {
                    break; // Not ready to print min_batch yet
                }
            }
        }
    }
}

fn lookup_parallel(n_threads: usize, query_path: &Path, index: SingleColoredKmers) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(2); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<QueryBatch>(2);

    let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
    let query_start = std::time::Instant::now();

    std::thread::scope(|s| {
        let reader_handle = s.spawn(|| {
            // Reader thread that pushes batches for workers

            let mut seq_id = 0_usize;
            while let Some(rec) = reader.read_next().unwrap() {
                // Todo: break long sequences into multiple batches and combine short sequences

                let mut seqs = SeqDB::new();
                seqs.push_seq(rec.seq);

                let batch = QueryBatch{
                    seqs,
                    seq_id_range: seq_id..seq_id+1,
                    start_kmer_offset: 0,
                    end_kmer_offset: kmers_in_n(index.k(), rec.seq.len()),
                    result: Vec::new(),
                };

                batch_send.send(batch).unwrap();

                seq_id += 1;
            }

            drop(batch_send); // Close the channel
        });

        let mut worker_handles = Vec::new();
        for _ in 0..n_threads {
            worker_handles.push(s.spawn(|| {
                while let Ok(mut batch) = batch_recv.recv() {
                    batch.run(&index);
                    output_send.send(batch).unwrap();
                }
            }));
        }

        // Wait for threads to finish
        reader_handle.join().unwrap(); // All work batches pushed to workers
        worker_handles.into_iter().for_each(|w| w.join().unwrap()); // All batches processed
    });

    let query_duration = query_start.elapsed();
    
    /*for color in 0..index.n_colors() {
        let hits = color_hit_counts[color];
        println!("Color {}: {} hits ({:.2}%)", color, hits, hits as f64 / total_kmers_queried as f64 * 100.0);
    }
    
    eprintln!("{} k-mers queried in {} seconds (excluding index loading time)", total_kmers_queried, query_duration.as_secs());
    eprintln!("{:.2}% of query k-mers found", color_hit_counts.iter().sum::<usize>() as f64 / total_kmers_queried as f64 * 100.0);
    eprintln!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / total_kmers_queried as f64);
    */

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

        Subcommands::Lookup{query: query_path, index: index_path, n_threads} => {
            eprintln!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(index_path).unwrap());

            let index_loading_start = std::time::Instant::now();
            let index = SingleColoredKmers::load(&mut index_input);
            eprintln!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());
            eprintln!("Running queries from {} ...", query_path.display());
            lookup_parallel(n_threads, &query_path, index);

        }
    } 
}