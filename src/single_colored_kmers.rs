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
use crate::lca_tree::LcaTree;
use crate::traits::*;

pub struct ColorStats {
    pub colored: usize,
    pub uncolored: usize,
    pub color_run_min: usize,
    pub color_run_max: usize,
    pub color_run_mean: f64,
}

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

#[derive(Debug, Clone)]
pub struct SingleColoredKmers<L: ContractLeft + Clone + MySerialize + From<LcsArray>, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: L,
    colors: C,
    color_names: Vec<String>,
    color_hierarchy: LcaTree,
}

impl<L: sbwt::ContractLeft + Clone + MySerialize + From<sbwt::LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize+ From<SimpleColorStorage>> ColoredKmerLookupAlgorithm for SingleColoredKmers<L, C> {
    fn lookup_kmers(&self, query: &[u8], k: usize) -> impl Iterator<Item = Option<usize>> {
        self.lookup_kmers(query, k)
    }
}

pub struct KmerLookupIterator<'a, 'b, L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, L>,
    index: &'a SingleColoredKmers<L, C>,
    query_pattern_length: usize,
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Iterator for KmerLookupIterator<'_, '_, L, C> {
    type Item = Option<usize>; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;
        let lcs = &self.index.lcs;
        let hierarchy = self.index.color_hierarchy();
        let root_id = hierarchy.root();

        if len >= self.query_pattern_length {
            // k-mer is found in the sbwt
            assert!(range.len() > 0);
            let mut color = self.index.get_color_of_range(range.clone());
            if color == Some(root_id) { return Some(color) }

            // Expand the interval as long as it has at most single color and the
            // LCS is at least the query pattern length.
            let (mut new_start, mut new_end) = (range.start, range.end);
            //eprintln!("Expanding from {}..{}, (len {})", new_start, new_end, new_end - new_start);
            while new_start > 0 && lcs.get_lcs(new_start) >= self.query_pattern_length {
                new_start -= 1;
                color = hierarchy.lca_options(color, self.index.get_color(new_start));
                if color == Some(root_id) { return Some(color) } // This is a Some(Some(color)). Means that the iterator produced something.
            }
            let n = self.index.sbwt.n_sets();
            while new_end < n && lcs.get_lcs(new_end) >= self.query_pattern_length {
                color = hierarchy.lca_options(color, self.index.get_color(new_end));
                if color == Some(root_id) { return Some(color) } // This is a Some(Some(color)). Means that the iterator produced something.
                new_end += 1;
            }

            Some(color) // This is a Some(Some(color)). Means that the iterator produced something.
        } else {
            Some(None) // Iterator not finished but the k-mer is not found -> no color
        }
    }
}


impl<T: AtomicUint> AtomicColorVec for Vec<T> {
    fn new(len: usize) -> Self {
        (0..len).map(|_| T::new(T::max_value())).collect()
    }

    fn update(&self, i: usize, x: usize, lca: &LcaTree) {
        assert!(x != Self::none_sentinel(), "x must not be the none sentinel");
        self[i].fetch_update(|cur| {
            if cur == Self::none_sentinel() {
                x               // first assignment: replace none with color
            } else {
                lca.lca(cur, x) // subsequent: merge via LCA
            }
        });
    }

    fn read(&self, i: usize) -> Option<usize> {
        let x = self[i].load(std::sync::atomic::Ordering::Relaxed);
        if x == Self::none_sentinel() {
            None
        } else {
            Some(x)
        }
    }

    fn none_sentinel() -> usize {
        T::max_value()
    }
}

struct ColoringBatch {
    dbs: Vec<(usize, SeqDB)>, // Pairs (color, seqs)
    dummy_mer_dbs: Vec<(usize, SeqDB)>, // Pairs (color, seqs)
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

    fn push_dummy_mer(&mut self, color: usize, mer: &[u8]) {
        // logic: push to the last DB if it has the right color,
        // othewise create a new DB and push to that. Add the length
        // of seq to self.total_len.

        let mut extended = false;
        if let Some((last_color, last_db)) = self.dummy_mer_dbs.last_mut() {
            if *last_color == color {
                // Extend last db
                last_db.push_seq(mer);
                extended = true;
            }
        }

        if !extended {
            // Start a new one
            let mut db = SeqDB::new();
            db.push_seq(mer);
            self.dummy_mer_dbs.push((color, db));
        }

        self.total_len += mer.len();
    }

    fn run<V: AtomicColorVec>(&self, si: &StreamingIndex<'_, SbwtIndex<SubsetMatrix>, LcsArray>, color_ids: &V, progress_counter: &AtomicU64, lca_tree: &LcaTree) {
        let k = si.k();
        let mut thread_progress = 0_usize;

        // Process full k-mers
        for (color, db) in self.dbs.iter() {
            for rec in db.iter() {
                let seq = rec.seq;
                let ms = si.matching_statistics_iter(seq);
                ms.enumerate().for_each(|(i, (len, range))| {
                    if len == k {
                        debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
                        color_ids.update(range.start, *color, lca_tree);
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

        // Process dummy mers
        for (color, db) in self.dummy_mer_dbs.iter() {
            for rec in db.iter() {
                let mer = rec.seq;
                let ms = si.matching_statistics_iter(mer);
                ms.for_each(|(len, range)| { 
                    assert!(range.len() > 0);
                    assert!(len < k);
                    // IMPORTANT. We assume that all the dummy paths to the starts of
                    // all the runs of ACGT are in the SBWT. Then this colex position
                    // will correspond to a dummy node. We could add a check for this
                    // but it's expensive.
                    let colex = range.start;
                    color_ids.update(colex, *color, lca_tree);
                });
            }
        }

        progress_counter.fetch_add(thread_progress as u64, std::sync::atomic::Ordering::Relaxed); // Add leftover progress
    }
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> SingleColoredKmers<L, C> {

    pub fn into_parts(self) -> (SbwtIndex<SubsetMatrix>, L, C, Vec<String>, LcaTree) {
        (self.sbwt, self.lcs, self.colors, self.color_names, self.color_hierarchy)
    }

    pub fn color_hierarchy(&self) -> &LcaTree {
        &self.color_hierarchy
    }

    pub fn color_names(&self) -> &[String] {
        &self.color_names
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
        5_u32
    }

    pub fn serialize(&self, mut out: &mut impl Write) {
        out.write_all(&Self::serialization_magic_constant()).unwrap();
        out.write_all(&Self::serialization_version_number().to_le_bytes()).unwrap();

        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out);

        self.colors.serialize(&mut out);
        self.color_hierarchy.serialize(&mut out).unwrap();

        for name in self.color_names.iter() {
            // Serialize length and bytes
            let name_bytes = name.as_bytes();
            bincode::serialize_into(&mut out, &(name_bytes.len() as u64)).unwrap();
            out.write_all(name_bytes).unwrap();
        }
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
        let color_hierarchy = LcaTree::load(input).unwrap();

        let mut color_names = Vec::new();
        for _ in 0..color_hierarchy.n_nodes() {
            let name_len: u64 = bincode::deserialize_from(&mut input).unwrap();
            let mut name_bytes = vec![0_u8; name_len as usize];
            input.read_exact(&mut name_bytes).unwrap();
            color_names.push(String::from_utf8(name_bytes).unwrap());
        }


        SingleColoredKmers{sbwt, lcs, colors, color_names, color_hierarchy}
    }

    pub fn n_colors_in_hierarchy(&self) -> usize {
        self.color_hierarchy.n_nodes()
    }

    pub fn n_kmers(&self) -> usize {
        self.sbwt.n_kmers()
    }

    fn color_run_stats(&self) -> (usize, usize, f64) { // (min, max, mean)
        let n = self.sbwt.n_sets();
        if n == 0 { return (0, 0, 0.0); }

        let mut min = usize::MAX;
        let mut max = 0_usize;
        let mut n_runs = 0_usize;
        let mut prev = self.get_color(0);
        let mut run_len = 1_usize;
        for i in 1..n {
            let cur = self.get_color(i);
            if cur == prev {
                run_len += 1;
            } else {
                min = min.min(run_len);
                max = max.max(run_len);
                n_runs += 1;
                run_len = 1;
                prev = cur;
            }
        }
        min = min.min(run_len);
        max = max.max(run_len);
        n_runs += 1;
        let mean = n as f64 / n_runs as f64;
        (min, max, mean)
    }

    pub fn color_stats(&self) -> ColorStats {
        let mut uncolored = 0_usize;
        let mut colored = 0_usize;
        for i in 0..self.sbwt.n_sets() {
            match self.get_color(i) {
                Some(_) => colored += 1,
                None => uncolored += 1,
            }
        }
        let (color_run_min, color_run_max, color_run_mean) = self.color_run_stats();
        ColorStats { colored, uncolored, color_run_min, color_run_max, color_run_mean }
    }

    pub fn get_color(&self, colex: usize) -> Option<usize> {
        assert!(colex <= self.sbwt.n_sets());
        self.colors.get_color(colex)
    }
    
    pub fn get_color_of_range(&self, colex_range: Range<usize>) -> Option<usize> {
        assert!(colex_range.end <= self.sbwt.n_sets());
        self.colors.get_color_of_range(colex_range, &self.color_hierarchy)
    }

    // Returns an iterator giving the color of each of the n-k+1 k-mers of the query.
    // k must be less or equal to the k in the SBWT index.
    // If query is shorter than k, returns an empty iterator.
    pub fn lookup_kmers<'a, 'b>(&'a self, query: &'b [u8], k: usize) -> KmerLookupIterator<'a, 'b, L, C>{
        assert!(k <= self.sbwt.k());
        let si = StreamingIndex {
            extend_right: &self.sbwt, 
            contract_left: &self.lcs,
            n: self.sbwt.n_sets(),
            k: self.sbwt.k(),
        };

        //let mut ms_iter = si.bounded_matching_statistics_iter(query, k);
        let mut ms_iter = si.matching_statistics_iter(query);

        // Skip over the first k-1 positions
        for _ in 0..k-1 {
            ms_iter.next(); // If the iterator ends early, will keep returning None
        }
        KmerLookupIterator { matching_stats_iter: ms_iter, index: self, query_pattern_length: k}

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
    fn mark_colors<T: SeqStream + Send, A: AtomicColorVec + Send + Sync>(sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray, input_streams: Vec<T>, n_threads: usize, color_hierarchy: &LcaTree) -> SimpleColorStorage {

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
                let mut batch = ColoringBatch{dbs: vec![], dummy_mer_dbs: vec![], total_len: 0};
                let b = 10000; // Batch size
                let k = sbwt.k();
                for (color, mut stream) in input_streams.into_iter().enumerate() {
                    while let Some(seq) = stream.stream_next() {
                        // Push dummy-mers
                        crate::util::for_each_run_with_key(seq, |c| IS_DNA[*c as usize], |mut run_range: Range<usize>| {
                            if !run_range.is_empty() && IS_DNA[seq[run_range.start] as usize] {
                                // Start of run of DNA characters
                                if run_range.len() >= k { // Clip to length k-1
                                    run_range = run_range.start..run_range.start+(k-1);
                                }
                                let mer = &seq[run_range.clone()];
                                batch.push_dummy_mer(color, mer);
                            }
                        });

                        // Push the rest in pieces, flushing when appropriate
                        crate::util::process_kmers_in_pieces(seq, k, b, |_piece_idx, piece: &[u8]|{
                            batch.push(color, piece);

                            if batch.total_len >= b {
                                // Swap the current batch with an empty batch, and send it to processing
                                let mut batch_to_send = ColoringBatch{dbs: vec![], dummy_mer_dbs: vec![], total_len: 0}; // Empty batch
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
                        batch.run(si_ref, color_ids_ref, n_bases_processed_ref, color_hierarchy);
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
        let mut total_some_count = 0_usize;
        let mut total_none_count = 0_usize;
        for i in 0..sbwt.n_sets() {
            let cv = color_ids.read(i);
            match cv {
                Some(_) => total_some_count += 1,
                None => total_none_count += 1,
            }
            compressed_colors.set_color(i, cv);
        }

        log::info!("Colored {} sbwt positions", total_some_count);
        log::info!("{} sbwt positions left uncolored", total_none_count);

        compressed_colors

    }


    pub fn new<T: SeqStream + Send>(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, input_streams: Vec<T>, mut color_names: Vec<String>, n_threads: usize, color_hierarchy: Option<LcaTree>) -> Self {
        assert!(color_names.len() == input_streams.len());
        for name in &color_names {
            assert!(name != "root" && name != "none", "color name '{}' is reserved", name);
        }

        // Build the color hierarchy before mark_colors so LCA-based merging can be used during indexing.
        let color_hierarchy = color_hierarchy.unwrap_or_else(|| {
            color_names.push("root".to_string());
            let n = color_names.len(); // n_colors + 1
            let edges = (0..n - 1).map(|i| (i, n - 1)).collect();
            LcaTree::new(n, edges).expect("default hierarchy cannot fail")
        });

        let required_bit_width = SimpleColorStorage::required_bit_width(color_hierarchy.n_nodes() + 1); // +1 for the "none"

        log::info!("Marking colors");
        let color_storage = if required_bit_width <= 8 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU8>>(&sbwt, &lcs, input_streams, n_threads, &color_hierarchy)
        } else if required_bit_width <= 16 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU16>>(&sbwt, &lcs, input_streams, n_threads, &color_hierarchy)
        } else if required_bit_width <= 32 {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU32>>(&sbwt, &lcs, input_streams, n_threads, &color_hierarchy)
        } else {
            SingleColoredKmers::<L,C>::mark_colors::<T, Vec::<AtomicU64>>(&sbwt, &lcs, input_streams, n_threads, &color_hierarchy)
        };

        Self::new_given_coloring(sbwt, lcs, color_storage, color_names, Some(color_hierarchy))
    }

    pub fn new_given_coloring(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, coloring: SimpleColorStorage, mut color_names: Vec<String>, color_hierarchy: Option<LcaTree>) -> Self {
        let n_colors = coloring.n_colors();
        for name in color_names.iter().take(n_colors) {
            assert!(name != "root" && name != "none", "color name '{}' is reserved", name);
        }

        let color_hierarchy = color_hierarchy.unwrap_or_else(|| {
            color_names.push("root".to_string());
            let n = color_names.len(); // n_colors + 1
            let edges = (0..n - 1).map(|i| (i, n - 1)).collect();
            LcaTree::new(n, edges).expect("default color hierarchy construction cannot fail")
        });

        log::info!("Indexing LCS array");
        let lcs_index = L::from(lcs);

        log::info!("Building Color id wavelet tree");
        let colors_index = C::from(coloring);

        log::info!("Color structure construction complete");
        SingleColoredKmers::<L, C> {
            sbwt, lcs: lcs_index, colors: colors_index, color_names, color_hierarchy,
        }
    }

    pub fn turn_nones_to_roots(&mut self) {
        let root_id = self.color_hierarchy.root();
        for i in 0..self.sbwt.n_sets() {
            if self.get_color(i) == None {
                self.colors.set_color(i, Some(root_id));
            }
        }
    }
}

pub struct SingleColoredKmersShort<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    inner: SingleColoredKmers<L,C>, // k-mers sharing an s-mer have been made to have the same color: the LCA in the color hierarchy
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> SingleColoredKmersShort<L,C> {

    // s is the query length. s <= k
    pub fn new(mut inner: SingleColoredKmers<L, C>, s: usize) -> Self {
        assert!(s <= inner.sbwt.k());
        for colex in 1..inner.sbwt.n_sets() {
            if inner.lcs.get_lcs(colex) >= s {
                let prev_color = inner.colors.get_color(colex-1);
                let now_color = inner.colors.get_color(colex);
                let updated = inner.color_hierarchy.lca_options(prev_color, now_color);
                inner.colors.set_color(colex-1, updated);
                inner.colors.set_color(colex, updated);
            }
        }
        Self { inner }
    }
}

impl<L: sbwt::ContractLeft + Clone + MySerialize + From<sbwt::LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize+ From<SimpleColorStorage>> ColoredKmerLookupAlgorithm for SingleColoredKmersShort<L, C> {
    fn lookup_kmers(&self, query: &[u8], s: usize) -> impl Iterator<Item = Option<usize>> {
        assert!(s <= self.inner.sbwt.k());
        let si = StreamingIndex {
            extend_right: &self.inner.sbwt, 
            contract_left: &self.inner.lcs,
            n: self.inner.sbwt.n_sets(),
            k: self.inner.sbwt.k(),
        };

        //let mut ms_iter = si.bounded_matching_statistics_iter(query, k);
        let mut ms_iter = si.matching_statistics_iter(query);

        // Skip over the first s-1 positions
        for _ in 0..s-1 {
            ms_iter.next(); // If the iterator ends early, will keep returning None
        }
        KmerLookupIteratorShort { matching_stats_iter: ms_iter, index: self, query_pattern_length: s}
    }
}

pub struct KmerLookupIteratorShort<'a, 'b, L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, L>,
    index: &'a SingleColoredKmersShort<L, C>,
    query_pattern_length: usize,
}

impl<L: ContractLeft + Clone + MySerialize + From<LcsArray> + LcsAccess, C: ColorStorage + Clone + MySerialize + From<SimpleColorStorage>> Iterator for KmerLookupIteratorShort<'_, '_, L, C> {
    type Item = Option<usize>; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len >= self.query_pattern_length {
            // s-mer is found in the sbwt
            assert!(range.len() > 0);

            // Because we have preprocessed to color array, all positions in the colex
            // range of an s-mer have the same color. So we can pick any of those. Let's
            // pick the first one.
            let color = self.index.inner.get_color(range.start);

            Some(color) // This is a Some(Some(color)). Means that the iterator produced something.
        } else {
            Some(None) // Iterator not finished but the k-mer is not found -> no color
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

impl LcsAccess for LcsWrapper {
    fn get_lcs(&self, colex: usize) -> usize {
        self.inner.access(colex)
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
