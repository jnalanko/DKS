#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::{collections::HashMap, fs::File, io::{BufRead, BufReader, BufWriter}, path::PathBuf};
use clap::{Parser, Subcommand};
use io::LazyFileSeqStream;
use jseqio::reader::DynamicFastXReader;
use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray};
use single_colored_kmers::SingleColoredKmers;
use parallel_queries::{TsvWriter, BedWriter};

mod single_colored_kmers;
mod io;
mod parallel_queries;
mod single_threaded_queries;
mod util;


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
        #[arg(short, required = true)]
        k: usize,

        #[arg(help = "A file with one fasta/fastq filename per line", short, long, required = true, help_heading = "Input")]
        input: PathBuf,

        #[arg(help = "Optional: a fasta/fastq file containing the unitigs of all the k-mers in the input files. More generally, any sequence file with same k-mers will do (unitigs, matchtigs, eulertigs...). This speeds up construction and reduces the RAM and disk usage", short, long, help_heading = "Input")]
        unitigs: Option<PathBuf>,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.", long = "external-memory")]
        temp_dir: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", short = 'f', long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Optional: a precomputed Bit Matrix SBWT file of the input k-mers.", short, long, help_heading = "Advanced use")]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Optional: a precomputed LCS file of the optional SBWT file.", short, long, help_heading = "Advanced use")]
        lcs_path: Option<PathBuf>,

    },

    #[command(arg_required_else_help = true)]
    Lookup {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Output in BED format instead of TSV. Requires --colors.", long, requires = "colors")]
        bed: bool,

        #[arg(help = "A tab-separated file with two columns: color_rank and color_name. Required when --bed is set.", long, requires = "bed")]
        colors: Option<PathBuf>,
    },

    #[command(arg_required_else_help = true, about = "Simple reference implementation for debugging this program.")]
    LookupDebug {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },
}

struct DynamicFastXReaderWrapper {
    inner: DynamicFastXReader,
}

impl sbwt::SeqStream for DynamicFastXReaderWrapper{
    fn stream_next(&mut self) -> Option<&[u8]> {
        self.inner.read_next().unwrap().map(|x| x.seq)
    }
}

fn load_bed_metadata(query_path: &PathBuf, colors_path: &PathBuf) -> (Vec<String>, HashMap<usize, String>) {

    // First pass: collect sequence names from query FASTA/FASTQ
    log::info!("Collecting sequence names from {} ...", query_path.display());
    let mut name_reader = DynamicFastXReader::from_file(&query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
    let mut seq_names = Vec::new();
    while let Some(rec) = name_reader.read_next().unwrap() {
        let header = std::str::from_utf8(rec.head).unwrap();
        // Take the first whitespace-delimited token as the sequence name
        let name = header.split_whitespace().next().unwrap_or(header);
        seq_names.push(name.to_string());
    }
    log::info!("Collected {} sequence names from query file", seq_names.len());

    // Parse colors file (tab-separated: color_rank\tcolor_name)
    let mut color_names = HashMap::new();
    let colors_file = BufReader::new(File::open(colors_path)
        .unwrap_or_else(|e| panic!("Could not open colors file {}: {e}", colors_path.display())));
    for line in colors_file.lines() {
        let line = line.unwrap();
        let line = line.trim();
        if line.is_empty() { continue; }
        let mut fields = line.split('\t');
        let rank: usize = fields.next()
            .unwrap_or_else(|| panic!("Missing color_rank in colors file"))
            .parse()
            .unwrap_or_else(|e| panic!("Invalid color_rank in colors file: {e}"));
        let name = fields.next()
            .unwrap_or_else(|| panic!("Missing color_name in colors file"))
            .to_string();
        color_names.insert(rank, name);
    }
    log::info!("Loaded {} color names from {}", color_names.len(), colors_path.display());

    (seq_names, color_names)
}

fn main() {

    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    log::info!("Running dks version {}", env!("CARGO_PKG_VERSION"));

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { input: input_fof, unitigs: unitigs_path, output: out_path, temp_dir, k, n_threads, forward_only, sbwt_path, lcs_path} => {
            let input_paths: Vec<PathBuf> = BufReader::new(File::open(&input_fof)
                .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", input_fof.display())))
                .lines().map(|f| PathBuf::from(f.unwrap())).collect();

            // Create output directory if does not exist
            std::fs::create_dir_all(out_path.parent().unwrap()).unwrap();
            let mut out = BufWriter::new(File::create(out_path.clone())
                .unwrap_or_else(|e| panic!("Could not create output file {}: {e}", out_path.display())));

            let add_rev_comps = !forward_only;

            let (sbwt, lcs) = if let Some(sbwt_path) = sbwt_path {
                let mut input = BufReader::new(File::open(&sbwt_path)
                    .unwrap_or_else(|e| panic!("Could not open SBWT file {}: {e}", sbwt_path.display())));
                let sbwt::SbwtIndexVariant::SubsetMatrix(sbwt) = sbwt::load_sbwt_index_variant(&mut input).unwrap();
                log::info!("Loaded SBWT with {} k-mers", sbwt.n_kmers());
                let lcs = if let Some(lcs_path) = lcs_path {
                    LcsArray::load(&mut BufReader::new(File::open(&lcs_path)
                        .unwrap_or_else(|e| panic!("Could not open LCS file {}: {e}", lcs_path.display())))).unwrap()
                } else {
                    LcsArray::from_sbwt(&sbwt, n_threads)
                };
                if sbwt.k() != k {
                    panic!("The k specified ({}) does not match the k of the provided SBWT ({})", k, sbwt.k());
                }
                (sbwt, lcs)
            } else {
                let all_input_seqs = if let Some(unitigs_path) = unitigs_path {
                    io::ChainedInputStream::new(vec![unitigs_path.clone()])
                } else {
                    io::ChainedInputStream::new(input_paths.clone())
                };
                let (sbwt, lcs) = if let Some(td) = temp_dir {
                    // Use disk-based construction
                    sbwt::SbwtIndexBuilder::new()
                        .add_rev_comp(add_rev_comps)
                        .k(k)
                        .build_lcs(true)
                        .n_threads(n_threads)
                        .precalc_length(8)
                        .algorithm(BitPackedKmerSortingDisk::new().dedup_batches(false).temp_dir(&td))
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
                (sbwt,lcs)
            };

            let individual_streams = input_paths.iter().map(|p| LazyFileSeqStream::new(p.clone(), add_rev_comps)).collect();
            log::info!("Marking colors");
            let index = SingleColoredKmers::new(sbwt, lcs, individual_streams, n_threads);

            log::info!("Writing to {}", out_path.display());
            index.serialize(&mut out);
            let out_size = std::fs::metadata(out_path).unwrap().len() as f64;
            log::info!("Index size on disk: {}" , human_bytes::human_bytes(out_size));
        },

        Subcommands::Lookup{query: query_path, index: index_path, n_threads, bed, colors} => {
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = SingleColoredKmers::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            // 128kb = 2^17 byte buffer
            let stdout = BufWriter::with_capacity(1 << 17, std::io::stdout());

            if bed {
                let colors_path = colors.unwrap(); // Safe: clap ensures --colors is present when --bed is set
                let (seq_names, color_names) = load_bed_metadata(&query_path, &colors_path);

                log::info!("Running queries from {} ...", query_path.display());
                let reader = DynamicFastXReader::from_file(&query_path)
                    .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
                let reader = DynamicFastXReaderWrapper { inner: reader };

                let writer = BedWriter::new(stdout, seq_names, color_names);
                parallel_queries::lookup_parallel(n_threads, reader, &index, 10000, writer);
            } else {
                log::info!("Running queries from {} ...", query_path.display());
                let reader = DynamicFastXReader::from_file(&query_path)
                    .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
                let reader = DynamicFastXReaderWrapper { inner: reader };

                let writer = TsvWriter::new(stdout);
                parallel_queries::lookup_parallel(n_threads, reader, &index, 10000, writer);
            };
        },

        Subcommands::LookupDebug{query: query_path, index: index_path} => {
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = SingleColoredKmers::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());
            log::info!("Running query debug implementation for {} ...", query_path.display());
            single_threaded_queries::lookup_single_threaded(&query_path, &index);

        }
    } 
}