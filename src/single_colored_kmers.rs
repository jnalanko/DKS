use std::io::{Read, Write};
use std::ops::Range;
use std::sync::atomic::{AtomicU16, AtomicU32, AtomicU64, AtomicU8};
use std::time::{Duration, Instant};

use bitvec::prelude::*;
use bitvec_sds::traits::RandomAccessU32;
//use bitvec_sds::wavelet_tree::{SelectSupportBoth, WaveletTree};
use crossbeam::channel::{Receiver, RecvTimeoutError};
use jseqio::seq_db::SeqDB;
use sbwt::{ContractLeft, LcsArray, MatchingStatisticsIterator, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use crate::color_storage::SimpleColorStorage;
use crate::traits::*;

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

#[derive(Debug, Clone)]
pub struct SingleColoredKmers<L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: L,
    colors: C,
    n_colors: usize,
}

impl<L: sbwt::ContractLeft + Clone + MySerialize + From<sbwt::LcsArray>, C: ColorStorage + Clone + MySerialize+ From<SimpleColorStorage>> ColoredKmerLookupAlgorithm for SingleColoredKmers<L, C> {
    fn lookup_kmers(&self, query: &[u8], k: usize) -> impl Iterator<Item = ColorVecValue> {
        self.lookup_kmers(query, k)
    }
}

pub struct KmerLookupIterator<'a, 'b, L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, L>,
    index: &'a SingleColoredKmers<L, C>,
    length_bound: usize,
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Iterator for KmerLookupIterator<'_, '_, L, C> {
    type Item = ColorVecValue; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len == self.length_bound {
            // k-mer is found in the sbwt
            debug_assert!(self.length_bound < self.index.k() || range.len() == 1); // Full k-mer should have a singleton range
            Some(self.index.get_color(range))
        } else {
            Some(ColorVecValue::None) // Iterator not finished but the k-mer is not found -> no color
        }
    }
}


impl<T: AtomicUint> AtomicColorVec for Vec<T> {
    fn update(&self, i: usize, x: usize) {
        assert!(x < T::max_value()-1); // MAX and MAX - 1 are reserved values
        self[i].fetch_update(|cur| {
            if cur == T::max_value() || cur == x { // No value stored yet, or storing the same value
                x // Update to x
            } else { // Different value already stored
                T::max_value()-1 // update to multiple
            }
        });
    }

    fn read(&self, i: usize) -> ColorVecValue {
        let x = self[i].load(std::sync::atomic::Ordering::Relaxed);
        if x == T::max_value(){
            ColorVecValue::None
        } else if x == T::max_value()-1 {
            ColorVecValue::Multiple
        } else {
            ColorVecValue::Single(x)
        }
    }

    fn new(len: usize) -> Self {
        (0..len).map(|_| T::new(T::max_value())).collect()
    }

}

struct ColoringBatch {
    dbs: Vec<(usize, SeqDB)>, // Pairs (color, seqs)
    total_len: usize,
}

impl ColoringBatch {
    fn push(&mut self, color: usize, seq: &[u8]) {
        // logic: push to the last DB if it has the right color,
        // othewise create a new DB and push to that. Add the length
        // of seq to self.total_len.

        let mut extended = false;
        if let Some((last_color, last_db)) = self.dbs.last_mut() {
            if *last_color == color {
                // Extend last db
                last_db.push_seq(seq);
                extended = true;
            }
        }

        if !extended {
            // Start a new one
            let mut db = SeqDB::new();
            db.push_seq(seq);
            self.dbs.push((color, db));
        }

        self.total_len += seq.len();
    }

    fn run<V: AtomicColorVec>(&self, si: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>, color_ids: &V, progress_counter: &AtomicU64) {
        let k = si.k();
        let mut thread_progress = 0_usize;
        for (color, db) in self.dbs.iter() {
            for rec in db.iter() {
                let seq = rec.seq;
                let ms = si.matching_statistics_iter(seq);
                ms.enumerate().for_each(|(i, (len, range))| { 
                    if len == k {
                        debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
                        color_ids.update(range.start, *color);
                    } else if cfg!(debug_assertions) && i >= k-1 {
                        // All valid k-mers should be found. If we're here, the k-mer must have had non-ACGT
                        // characters which make it invalid. Let's verify that.
                        let kmer = &seq[i-(k-1)..=i];
                        let all_ACGT = kmer.iter().all(|c| IS_DNA[*c as usize]);
                        if all_ACGT {
                            panic!("Error: k-mer {} not found in sbwt", String::from_utf8_lossy(kmer));
                        }
                    }
                    thread_progress += 1;
                    if thread_progress == 10000 {
                        // Only record progress every 10000 iterations to reduce synchronization overhead. 
                        // This made the code 30% faster in tests with 4 threads.
                        progress_counter.fetch_add(10000, std::sync::atomic::Ordering::Relaxed);
                        thread_progress = 0;
                    }
                });
            }
        }

        progress_counter.fetch_add(10000, std::sync::atomic::Ordering::Relaxed);
    }
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> SingleColoredKmers<L, C> {

    pub fn into_parts(self) -> (SbwtIndex<SubsetMatrix>, L, C, usize) {
        (self.sbwt, self.lcs, self.colors, self.n_colors)
    }

    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    // This is used to identify files serialized from this struct
    fn serialization_magic_constant() -> [u8; 4] {
        [54, 229, 250, 84] // Four randomly generated bytes
    }

    // This is used to identify the version of the serialization format
    fn serialization_version_number() -> u32 {
        2_u32
    }

    pub fn serialize(&self, mut out: &mut impl Write) {
        out.write_all(&Self::serialization_magic_constant()).unwrap();
        out.write_all(&Self::serialization_version_number().to_le_bytes()).unwrap();

        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out);

        self.colors.serialize(&mut out);
        bincode::serialize_into(&mut out, &self.n_colors).unwrap();
    }

    pub fn load(mut input: &mut impl Read) -> SingleColoredKmers<L, C> {

        let mut magic = [0_u8; 4];
        input.read_exact(&mut magic).unwrap();
        if magic != Self::serialization_magic_constant() {
            panic!("Error loading index: invalid file format (magic constant mismatch)");
        }

        let mut version_bytes = [0_u8; 4];
        input.read_exact(&mut version_bytes).unwrap();
        let version = u32::from_le_bytes(version_bytes);
        if version != Self::serialization_version_number() {
            panic!("Error loading index: wrong file format version number (found {}, expected {})", version, Self::serialization_version_number());
        }

        let sbwt = SbwtIndex::<sbwt::SubsetMatrix>::load(input).unwrap();
        let lcs = *L::load(input);

        let colors = *C::load(&mut input);
        let n_colors = bincode::deserialize_from(&mut input).unwrap();

        SingleColoredKmers{sbwt, lcs, colors, n_colors}
    }

    #[allow(dead_code)]
    pub fn n_colors(&self) -> usize {
        self.n_colors
    }
    
    pub fn get_color(&self, colex_range: Range<usize>) -> ColorVecValue {
        assert!(colex_range.end <= self.sbwt.n_sets());
        self.colors.get_color_of_range(colex_range)
    }

    // Returns an iterator giving the color of each of the n-k+1 k-mers of the query.
    // k must be less or equal to the k in the SBWT index.
    // If query is shorter than k, returns an empty.
    pub fn lookup_kmers<'a, 'b>(&'a self, query: &'b [u8], k: usize) -> KmerLookupIterator<'a, 'b, L, C>{
        assert!(k <= self.sbwt.k());
        let si = StreamingIndex {
            extend_right: &self.sbwt, 
            contract_left: &self.lcs,
            n: self.sbwt.n_sets(),
            k: self.sbwt.k(), // Can be different than the k given as a parameter to this function!
        };

        let mut ms_iter = si.bounded_matching_statistics_iter(query, k);

        // Skip over the first k-1 positions
        for _ in 0..k-1 {
            ms_iter.next(); // If the iterator ends early, will keep returning None
        }
        KmerLookupIterator { matching_stats_iter: ms_iter, index: self, length_bound: k}

    }

    fn progress_print_thread(n_bases_processed: &AtomicU64, quit_signal: Receiver<bool>) {
        log::info!("Processing up to 2n bases, where n is the number of bases in the input."); // 2n due to reverse complements
        let print_interval = 10; // seconds
        let mut last_print_time = Instant::now();
        let mut last_count = 0_u64;
        loop {
            match quit_signal.recv_timeout(Duration::from_secs(print_interval)) {
                Ok(_) => return, // Received quit signal
                Err(RecvTimeoutError::Timeout) => { // print_interval seconds has passed
                    let count = n_bases_processed.load(std::sync::atomic::Ordering::Relaxed);
                    let dcount = count - last_count;
                    let elapsed = last_print_time.elapsed().as_secs_f64();
                    log::info!("{} bases processed (total {}) ({:.2} Mbases/sec)", dcount, count, dcount as f64 / elapsed / 1e6);
                    last_count = count;
                    last_print_time = Instant::now();
                },
                Err(RecvTimeoutError::Disconnected) => {
                    // I'm not sure when this would happen, but let's just quit
                    return 
                }
            }
        }
    }


    // Generic function that works on any of u8, 16, u32 and u64
    fn mark_colors<T: SeqStream + Send, A: AtomicColorVec + Send + Sync>(sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray,  input_streams: Vec<T>, n_threads: usize) -> SimpleColorStorage {

        let color_ids = A::new(sbwt.n_sets());
        let si = StreamingIndex::new(sbwt, lcs);
        let n_colors = input_streams.len();

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        let n_bases_processed = AtomicU64::new(0);
        std::thread::scope(|scope| { thread_pool.install(|| {

            // Let's set up a thread that prints progress in regular intervals.
            // This channel will be used to tell it to quit:
            let (quit_print_send, quit_print_recv) = crossbeam::channel::unbounded::<bool>();
            let _progress_printer = scope.spawn({
                let n_bases_processed = &n_bases_processed;
                move || {
                    SingleColoredKmers::<L,C>::progress_print_thread(n_bases_processed, quit_print_recv);
                }
            });

            let (batch_send, batch_recv) = crossbeam::channel::bounded::<ColoringBatch>(4);

            // Create a reader that pushes batches to workers
            let reader_handle = scope.spawn(move || {
                let mut batch = ColoringBatch{dbs: vec![], total_len: 0};
                let b = 10000; // Batch size
                let k = sbwt.k();
                for (color, mut stream) in input_streams.into_iter().enumerate() {
                    while let Some(seq) = stream.stream_next() {
                        crate::util::process_kmers_in_pieces(seq, k, b, |_piece_idx, piece: &[u8]|{
                            batch.push(color, piece);

                            if batch.total_len >= b {
                                // Swap the current batch with an empty batch, and send it to processing
                                let mut batch_to_send = ColoringBatch{dbs: vec![], total_len: 0}; // Empty batch
                                std::mem::swap(&mut batch, &mut batch_to_send);
                                batch_send.send(batch_to_send).unwrap();
                            }
                        });
                    }
                }
                // Push the last batch
                if batch.total_len > 0 {
                    batch_send.send(batch).unwrap();
                }

                // batch_send is dropped here which closes the channel
            });

            // Create worker threads
            let mut worker_handles = Vec::new();
            for _ in 0..n_threads {
                let batch_recv_clone = batch_recv.clone(); // Moved into worker
                let si_ref = &si; // Moved into worker
                let n_bases_processed_ref = &n_bases_processed; // Moved into worker
                let color_ids_ref = &color_ids; // Moved into worker
                worker_handles.push(scope.spawn(move || {
                    while let Ok(batch) = batch_recv_clone.recv() {
                        batch.run(si_ref, color_ids_ref, n_bases_processed_ref);
                    }
                }));
            }

            // Wait for reader to finish
            reader_handle.join().unwrap();

            // Wait for the workers to finish.
            for w in worker_handles {
                w.join().unwrap();
            }

            // Tell the progress printer to quit (otherwise we hang)
            quit_print_send.send(true).unwrap(); 
        })});

        // Compress color_ids into a BitVec
        log::info!("Bitpacking color id array");
        let mut compressed_colors = SimpleColorStorage::new(sbwt.n_sets(), n_colors);
        let mut total_single_count = 0_usize;
        let mut total_multiple_count = 0_usize;
        for i in 0..sbwt.n_sets() {
            let color = color_ids.read(i);
            match color {
                ColorVecValue::Single(_) => total_single_count += 1,
                ColorVecValue::Multiple => total_multiple_count += 1,
                ColorVecValue::None => {},
            }
            compressed_colors.set_color(i, color); 

        }
        assert!(total_single_count + total_multiple_count <= sbwt.n_kmers());

        log::info!("Marked {} single-colored k-mers", total_single_count);
        log::info!("Marked {} multi-colored k-mers", total_multiple_count);
        log::info!("{} k-mers left not marked", sbwt.n_kmers() - total_single_count - total_multiple_count);

        compressed_colors

    }


    pub fn new<T: SeqStream + Send>(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, input_streams: Vec<T>, n_threads: usize) -> Self {

        let n_colors = input_streams.len();

        let required_bit_width = SimpleColorStorage::required_bit_width(n_colors);

        log::info!("Marking colors");
        let color_storage = if required_bit_width <= 8 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU8>>(&sbwt, &lcs, input_streams, n_threads)
        } else if required_bit_width <= 16 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU16>>(&sbwt, &lcs, input_streams, n_threads)
        } else if required_bit_width <= 32 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU32>>(&sbwt, &lcs, input_streams, n_threads)
        } else {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU64>>(&sbwt, &lcs, input_streams, n_threads)
        };

        Self::new_given_coloring(sbwt, lcs, color_storage)
    }

    pub fn new_given_coloring(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, coloring: SimpleColorStorage) -> Self {
        let n_colors = coloring.n_colors();

        log::info!("Indexing LCS array");
        let lcs_index = L::from(lcs);

        log::info!("Building Color id wavelet tree");
        let colors_index = C::from(coloring);

        log::info!("Color structure construction complete");
        SingleColoredKmers::<L, C> {
            sbwt, lcs: lcs_index, n_colors, colors: colors_index,
        }

    }
}

// Wrapper so that we can implement the foreing trait RandomAccessU32
#[derive(Debug, Clone)]
pub struct LcsWrapper {
    pub inner: LcsArray
}

impl ContractLeft for LcsWrapper {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        self.inner.contract_left(I, target_len)
    }
}

impl From<LcsArray> for LcsWrapper {
    fn from(lcs: LcsArray) -> Self {
        Self {inner: lcs}
    }
}

impl RandomAccessU32 for LcsWrapper {
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn get(&self, idx: usize) -> u32 {
        self.inner.access(idx) as u32
    }
}

impl MySerialize for LcsWrapper {
    fn serialize(&self, out: &mut impl Write) {
        self.inner.serialize(out).unwrap();
    }

    fn load(input: &mut impl Read) -> Box<Self> {
        let inner = LcsArray::load(input).unwrap();
        Box::new(LcsWrapper { inner })
    }
}
