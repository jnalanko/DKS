#![allow(non_snake_case, clippy::needless_range_loop, clippy::len_zero)] // Using upper-case variable names from the source material

use std::{fs::File, io::{BufRead, BufReader, BufWriter, Read, Write}, path::PathBuf};
use clap::{Parser, Subcommand};
use io::LazyFileSeqStream;
use jseqio::reader::DynamicFastXReader;
use sbwt::{BitPackedKmerSortingDisk, BitPackedKmerSortingMem, LcsArray};
use single_colored_kmers::SingleColoredKmers;
use parallel_queries::OutputWriter;

use crate::{color_storage::SimpleColorStorage, parallel_queries::RunWriter, single_colored_kmers::{ColorStats, LcsWrapper}, traits::ColorVecValue, wavelet_tree::WaveletTreeWrapper};

mod single_colored_kmers;
mod io;
mod parallel_queries;
mod single_threaded_queries;
mod util;
mod wavelet_tree;
mod traits;
mod color_storage;

type FixedKColorIndex = SingleColoredKmers<LcsWrapper, SimpleColorStorage>;
type FlexibleKColorIndex = SingleColoredKmers<WaveletTreeWrapper, SimpleColorStorage>;

enum ColorIndex {
    FixedK(FixedKColorIndex),
    FlexibleK(FlexibleKColorIndex),
}

const DKS_FILE_ID: [u8; 8] = *b"dks0.1.2";
const FIXED_INDEX_TYPE_ID: [u8; 4] = *b"fixd";
const FLEXIBLE_INDEX_TYPE_ID: [u8; 4] = *b"flex";

impl ColorIndex {
    fn serialize(&self, out: &mut impl Write) {
        match self {
            ColorIndex::FixedK(index) => {
                out.write_all(&DKS_FILE_ID).unwrap();
                out.write_all(&FIXED_INDEX_TYPE_ID).unwrap();
                index.serialize(out);
            },
            ColorIndex::FlexibleK(index) => {
                out.write_all(&DKS_FILE_ID).unwrap();
                out.write_all(&FLEXIBLE_INDEX_TYPE_ID).unwrap();
                index.serialize(out);
            }
        }
    }

    fn load(input: &mut impl Read) -> Self {
        let mut file_id = [0_u8; 8];
        input.read_exact(&mut file_id).unwrap();
        assert_eq!(file_id, DKS_FILE_ID, "Invalid DKS file ID");

        let mut type_id = [0_u8; 4];
        input.read_exact(&mut type_id).unwrap();
        match type_id {
            FIXED_INDEX_TYPE_ID => {
                let index = ColorIndex::FixedK(FixedKColorIndex::load(input));
                log::info!("Loaded fixed-k index with k = {}", index.k());
                index
            },
            FLEXIBLE_INDEX_TYPE_ID => {
                let index = ColorIndex::FlexibleK(FlexibleKColorIndex::load(input));
                log::info!("Loaded flexible-k index with k = {}", index.k());
                index
            },
            _ => {
                panic!("Unknown index type ID in DKS file: {}", String::from_utf8_lossy(&type_id));
            }
        }
    }

    fn k(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.k(),
            ColorIndex::FlexibleK(index) => index.k(),
        }
    }

    fn is_flexible(&self) -> bool {
        match self {
            ColorIndex::FlexibleK(_) => true,
            ColorIndex::FixedK(_) => false,
        }
    }

    fn color_names(&self) -> &[String] {
        match self {
            ColorIndex::FixedK(index) => index.color_names(),
            ColorIndex::FlexibleK(index) => index.color_names(),
        }
    }

    fn n_colors(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.n_colors(),
            ColorIndex::FlexibleK(index) => index.n_colors(),
        }
    }

    fn n_kmers(&self) -> usize {
        match self {
            ColorIndex::FixedK(index) => index.n_kmers(),
            ColorIndex::FlexibleK(index) => index.n_kmers(),
        }
    }

    fn color_stats(&self) -> ColorStats {
        match self {
            ColorIndex::FixedK(index) => index.color_stats(),
            ColorIndex::FlexibleK(index) => index.color_stats(),
        }
    }
}

fn into_flexible_index(fixed_index: FixedKColorIndex) -> FlexibleKColorIndex {
    let (sbwt, lcs, coloring, _n_colors, color_names) = fixed_index.into_parts();
    FlexibleKColorIndex::new_given_coloring(sbwt, lcs.inner, coloring, color_names)
}

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

        #[arg(help = "A file with one fasta/fastq filename per line, one per color", short, long, help_heading = "Input", conflicts_with = "sequence_colors")]
        file_colors: Option<PathBuf>,

        #[arg(help = "Give input as a single file, one sequence per color", short, long, help_heading = "Input", conflicts_with = "file_colors")]
        sequence_colors: Option<PathBuf>,

        #[arg(help = "Optional: a fasta/fastq file containing the unitigs of all the k-mers in the input files. More generally, any sequence file with same k-mers will do (unitigs, matchtigs, eulertigs...). This speeds up construction and reduces the RAM and disk usage", short, long, help_heading = "Input")]
        unitigs: Option<PathBuf>,

        #[arg(help = "Output filename", short, long, required = true)]
        output: PathBuf,

        #[arg(help = "Run in external memory construction mode using the given directory as temporary working space. This reduces the RAM peak but is slower. The resulting index will still be exactly the same.", long = "external-memory")]
        temp_dir: Option<PathBuf>,

        #[arg(help = "Do not add reverse complemented k-mers", long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Optional: a precomputed Bit Matrix SBWT file of the input k-mers. Must have been built with --add-all-dummy-paths", short = 'b', long, help_heading = "Advanced use")]
        sbwt_path: Option<PathBuf>,

        #[arg(help = "Optional: a precomputed LCS file of the optional SBWT file. Must have been built with --add-all-dummy-paths", short, long, help_heading = "Advanced use")]
        lcs_path: Option<PathBuf>,

        #[arg(help = "Build a flexible index supporting queries for any s-mer with s <= k. The index is slightly larger and the queries are approximately 3-10x slower.", long = "flexible", default_value = "false")]
        flexible: bool,

        #[arg(help = "Optional: a file with one color name per line, in the same order as the input files. Defaults to using the input filenames as color names.", long = "color-names", help_heading = "Input")]
        color_names_file: Option<PathBuf>,

    },

    #[command(arg_required_else_help = true)]
    Lookup {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,

        #[arg(help = "Query k-mer length. Must be less or equal to the k used in index construction. If not given, defaults to the same k as during index construction.", short, required = false)]
        k: Option<usize>,

        #[arg(help = "Print color names instead of color rank integers. K-mers present in multiple colors are reported as '*' normally, or 'multiple' when this flag is set.", long = "report-color-names")]
        report_color_names: bool,

        #[arg(help = "Print query names instead of query rank integers.", long = "report-query-names")]
        report_query_names: bool,

        #[arg(help = "Print lines for runs of k-mers not found in the index. The miss symbol is '-' normally, or 'none' when --report-color-names is set.", long = "report-misses")]
        report_misses: bool,

        #[arg(help = "Do not print the header line.", long = "no-header")]
        no_header: bool,
    },

    #[command(about = "Print statistics about an index file.")]
    Stats {
        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    #[command(arg_required_else_help = true, about = "Simple reference implementation for debugging this program.")]
    LookupDebug {
        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "Path to the index file", short, long, required = true)]
        index: PathBuf,
    },

    #[command(arg_required_else_help = true, about = "Debug: build individual SBWTs per color and query them separately.")]
    IndividualSbwtDebug {
        #[arg(help = "A file with one fasta/fastq filename per line (one per color)", short, long, required = true, help_heading = "Input")]
        input: PathBuf,

        #[arg(help = "A fasta/fastq query file", short, long, required = true)]
        query: PathBuf,

        #[arg(help = "K-mer length (used for both indexing and querying)", short, required = true)]
        k: usize,

        #[arg(help = "Do not add reverse complemented k-mers", short = 'f', long = "forward-only")]
        forward_only: bool,

        #[arg(help = "Number of parallel threads", short = 't', long = "n-threads", default_value = "4")]
        n_threads: usize,
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

fn load_seq_names(query_path: &PathBuf) -> Vec<String> {
    log::info!("Collecting sequence names from {} ...", query_path.display());
    let mut name_reader = DynamicFastXReader::from_file(query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));
    let mut seq_names = Vec::new();
    while let Some(rec) = name_reader.read_next().unwrap() {
        let header = std::str::from_utf8(rec.head).unwrap();
        let name = header.split_whitespace().next().unwrap_or(header);
        seq_names.push(name.to_string());
    }
    log::info!("Collected {} sequence names from query file", seq_names.len());
    seq_names
}

fn run_queries<W: RunWriter>(n_threads: usize, reader: DynamicFastXReader, index: ColorIndex, batch_size: usize, k: usize, writer: W) {
    let reader = DynamicFastXReaderWrapper { inner: reader };
    match index {
        ColorIndex::FixedK(index) => {
            if k < index.k() {
                log::warn!("Querying with k shorter than the k that was used for indexing ({} < {}). This will give the right answers, but might be very slow! Consider indexing with --flexible instead for better performance.", k, index.k());
            }
            parallel_queries::lookup_parallel(n_threads, reader, &index, batch_size, k, writer);
        },
        ColorIndex::FlexibleK(index) => {
            parallel_queries::lookup_parallel(n_threads, reader, &index, batch_size, k, writer);
        }
    }
}


fn individual_sbwt_debug(input_fof: &PathBuf, query_path: &PathBuf, k: usize, forward_only: bool, n_threads: usize) {
    // Read input file-of-files
    let input_paths: Vec<PathBuf> = BufReader::new(File::open(input_fof)
        .unwrap_or_else(|e| panic!("Could not open input file {}: {e}", input_fof.display())))
        .lines().map(|f| PathBuf::from(f.unwrap())).collect();

    // Build one SBWT per input file, keeping all in memory
    let mut indices = Vec::new();
    for input_path in &input_paths {
        log::info!("Building SBWT for {}", input_path.display());
        let stream = LazyFileSeqStream::new(input_path.clone(), false);
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(!forward_only)
            .k(k)
            .build_lcs(true)
            .n_threads(n_threads)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
            .run(stream);
        indices.push((sbwt, lcs.unwrap()));
    }

    // Open query file
    let mut reader = DynamicFastXReader::from_file(query_path)
        .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));

    // Create OutputWriter (same format as lookup command)
    let stdout = BufWriter::with_capacity(1 << 17, std::io::stdout());
    let mut writer = OutputWriter::new(stdout, None, None, false, true);
    writer.write_header();

    // Stream query sequences, querying all SBWTs in lockstep
    let mut seq_id: isize = 0;
    while let Some(rec) = reader.read_next().unwrap() {
        let mut ms_iters: Vec<_> = indices.iter()
            .map(|(sbwt, lcs)| {
                let si = sbwt::StreamingIndex::new(sbwt, lcs);
                si.matching_statistics_iter(rec.seq)
            })
            .collect();

        // Skip first k-1 positions (not full k-mers yet)
        for _ in 0..k.saturating_sub(1) {
            for iter in ms_iters.iter_mut() { iter.next(); }
        }

        // Advance all iterators in lockstep, run-length encoding on the fly
        let mut run_start = 0usize;
        let mut run_color = ColorVecValue::None;
        let mut kmer_count = 0usize;
        loop {
            let steps: Vec<Option<_>> = ms_iters.iter_mut().map(|it| it.next()).collect();
            if steps.iter().any(|s| s.is_none()) { break; }

            // Union colors: a k-mer is in color i iff its MS length == k
            let color = steps.into_iter().enumerate().fold(
                ColorVecValue::None,
                |acc, (color_id, step)| {
                    let (len, _) = step.unwrap();
                    if len == k { acc.union(ColorVecValue::Single(color_id)) }
                    else { acc }
                },
            );

            if kmer_count == 0 {
                run_start = 0;
                run_color = color;
            } else if color != run_color {
                writer.write_run(seq_id, run_color, run_start..kmer_count);
                run_start = kmer_count;
                run_color = color;
            }
            kmer_count += 1;
        }
        if kmer_count > 0 {
            writer.write_run(seq_id, run_color, run_start..kmer_count);
        }
        seq_id += 1;
    }
    writer.flush();
}

fn main() {

    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    log::info!("Running dks version {}", env!("CARGO_PKG_VERSION"));

    let args = Cli::parse();

    match args.command {
        Subcommands::Build { file_colors, sequence_colors, unitigs: unitigs_path, output: out_path, temp_dir, k, n_threads, forward_only, sbwt_path, lcs_path, flexible, color_names_file} => {
            let input_fof = if let Some(file_colors) = file_colors {
                file_colors
            } else {
                unimplemented!() // Sequence colors
            };
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
                        .add_all_dummy_paths(true) // This is required for multi-k support
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
                        .add_all_dummy_paths(true) // This is required for multi-k support
                        .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
                    .run(all_input_seqs)
                };
                let lcs = lcs.unwrap(); // Ok because of build_lcs(true)
                (sbwt,lcs)
            };

            let individual_streams = input_paths.iter().map(|p| LazyFileSeqStream::new(p.clone(), add_rev_comps)).collect();
            let color_names: Vec<String> = if let Some(ref names_path) = color_names_file {
                let names: Vec<String> = BufReader::new(File::open(names_path)
                    .unwrap_or_else(|e| panic!("Could not open color names file {}: {e}", names_path.display())))
                    .lines().map(|l| l.unwrap()).collect();
                if names.len() != input_paths.len() {
                    panic!("Color names file has {} names but there are {} input files", names.len(), input_paths.len());
                }
                names
            } else {
                input_paths.iter().map(|p| p.as_os_str().to_str().unwrap().to_owned()).collect()
            };
            log::info!("Marking colors");
            let index = FixedKColorIndex::new(sbwt, lcs, individual_streams, color_names, n_threads);

            let index = if flexible {
                log::info!("Transforming index to support flexible queries");
                ColorIndex::FlexibleK(into_flexible_index(index))
            } else {
                ColorIndex::FixedK(index)
            };

            log::info!("Writing to {}", out_path.display());
            index.serialize(&mut out);
            let out_size = std::fs::metadata(out_path).unwrap().len() as f64;
            log::info!("Index size on disk: {}" , human_bytes::human_bytes(out_size));

        },

        Subcommands::Lookup{query: query_path, index: index_path, n_threads, k, report_color_names, report_query_names, report_misses, no_header} => {
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = ColorIndex::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());

            let k = k.unwrap_or(index.k());
            if k > index.k() {
                panic!("Error: query k = {} larger than indexing k = {}", k, index.k());
            } else if k == index.k() && index.is_flexible() {
                log::warn!("Running with query k equal to indexing k. For faster queries, build a fixed-k index instead (no --flexible at indexing)");
            }

            let color_names = report_color_names.then(|| index.color_names().to_vec());
            let seq_names = report_query_names.then(|| load_seq_names(&query_path));

            let reader = DynamicFastXReader::from_file(&query_path)
                .unwrap_or_else(|e| panic!("Could not open query file {}: {e}", query_path.display()));

            let stdout = BufWriter::with_capacity(1 << 17, std::io::stdout());
            let writer = OutputWriter::new(stdout, seq_names, color_names, report_misses, !no_header);

            let batch_size = 10000;
            log::info!("Running queries from {} ...", query_path.display());
            run_queries(n_threads, reader, index, batch_size, k, writer);
        },

        Subcommands::Stats { index: index_path } => {
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));
            let index = ColorIndex::load(&mut index_input);

            let stats = index.color_stats();
            println!("Index type:            {}", if index.is_flexible() { "flexible-k" } else { "fixed-k" });
            println!("k:                     {}", index.k());
            println!("Number of colors:      {}", index.n_colors());
            println!("Number of k-mers:      {}", index.n_kmers());
            println!("Single-colored k-mers: {}", stats.single);
            println!("Multi-colored k-mers:  {}", stats.multiple);
            println!("Uncolored k-mers:      {}", stats.uncolored);
            println!("Color run min length:  {}", stats.color_run_min);
            println!("Color run max length:  {}", stats.color_run_max);
            println!("Color run mean length: {:.2}", stats.color_run_mean);
        },

        Subcommands::IndividualSbwtDebug { input: input_fof, query: query_path, k, forward_only, n_threads } => {
            individual_sbwt_debug(&input_fof, &query_path, k, forward_only, n_threads);
        },

        Subcommands::LookupDebug{query: query_path, index: index_path} => {
            log::info!("Loading the index ...");
            let mut index_input = BufReader::new(File::open(&index_path)
                .unwrap_or_else(|e| panic!("Could not open index file {}: {e}", index_path.display())));

            let index_loading_start = std::time::Instant::now();
            let index = ColorIndex::load(&mut index_input);
            log::info!("Index loaded in {} seconds", index_loading_start.elapsed().as_secs_f64());
            log::info!("Running query debug implementation for {} ...", query_path.display());

            match index {
                ColorIndex::FixedK(index) => {
                    single_threaded_queries::lookup_single_threaded(&query_path, &index, index.k());
                },
                ColorIndex::FlexibleK(index) => {
                    single_threaded_queries::lookup_single_threaded(&query_path, &index, index.k());
                }
            }

        }
    } 
}