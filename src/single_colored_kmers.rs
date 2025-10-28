use std::io::{Read, Write};
use std::sync::atomic::{AtomicU16, AtomicU32, AtomicU64, AtomicU8};
use std::time::{Duration, Instant};
use std::sync::atomic::Ordering::Relaxed;

use bitvec::prelude::*;
use bitvec::{field::BitField, order::Lsb0, vec::BitVec};
use crossbeam::channel::{Receiver, RecvTimeoutError};
use sbwt::{LcsArray, MatchingStatisticsIterator, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

#[derive(PartialEq, Eq, Debug, Copy, Clone)]
pub enum ColorVecValue {
    Single(usize),
    Multiple,
    None,
} 

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColorStorage {
    colors: BitVec::<usize, Lsb0>,
    bits_per_color: usize,
}

impl ColorStorage {
    fn get_color(&self, colex: usize) -> ColorVecValue {
        let x: usize = self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].load_le();
        if x == (1 << self.bits_per_color) - 1 { // Max value is reserved for None
            ColorVecValue::None
        } else if x == (1 << self.bits_per_color) - 2 { // Max - 1 is reserved for Multiple
            ColorVecValue::Multiple
        } else {
            ColorVecValue::Single(x)
        }
    }

    fn set_color(&mut self, colex: usize, color: ColorVecValue) {
        let value = match color {
            ColorVecValue::Single(x) => {
                assert!(x < (1 << self.bits_per_color) - 2); 
                x
            },
            ColorVecValue::Multiple => (1 << self.bits_per_color) - 2,
            ColorVecValue::None => (1 << self.bits_per_color) - 1,
        };
        self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].store_le(value);
    }

    fn new(len: usize, n_colors: usize) -> Self {
        let bits_per_color = Self::required_bit_width(n_colors);
        ColorStorage {
            colors: bitvec![0; len * bits_per_color],
            bits_per_color,
        }
    }

    fn required_bit_width(n_colors: usize) -> usize {
        log2_ceil(n_colors + 2) // +2 to reserve two special values: one for "no color" and one for "multiple colors"
    }
}

#[derive(Debug, Clone)]
pub struct SingleColoredKmers {
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    colors: ColorStorage,
    n_colors: usize,
}

pub struct KmerLookupIterator<'a, 'b> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, LcsArray>,
    index: &'a SingleColoredKmers,
}

impl Iterator for KmerLookupIterator<'_, '_> {
    type Item = ColorVecValue; // Color id of k-mer

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len == self.index.sbwt.k() {
            // k-mer is found in the sbwt
            debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
            let colex = range.start;
            Some(self.index.get_color(colex))
        } else {
            Some(ColorVecValue::None) // Iterator not finished but the k-mer is not found -> no color
        }
    }
}

trait AtomicColorVec{
    // Represents None as the max value of the atomic type

    fn update(&self, i: usize, x: usize); 
    fn read(&self, i: usize) -> ColorVecValue;
    fn new(len: usize) -> Self; // Stores a None (=max_value()) to each position
}

trait AtomicUint {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize;
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F);
    fn max_value() -> usize;
    fn new(val: usize) -> Self;
}

impl AtomicUint for AtomicU8 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u8)).unwrap();
    }
    fn max_value() -> usize {
        u8::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU8::new(val as u8)
    }
}

impl AtomicUint for AtomicU16 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u16)).unwrap();
    }
    fn max_value() -> usize {
        u16::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU16::new(val as u16)
    }
}

impl AtomicUint for AtomicU32 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u32)).unwrap();
    }
    fn max_value() -> usize {
        u32::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU32::new(val as u32)
    }
}

impl AtomicUint for AtomicU64 {
    fn load(&self, order: std::sync::atomic::Ordering) -> usize {
        self.load(order) as usize
    }
    fn fetch_update<F: Fn(usize) -> usize>(&self, f: F) {
        self.fetch_update(Relaxed, Relaxed, |x| Some(f(x as usize) as u64)).unwrap();
    }
    fn max_value() -> usize {
        u64::MAX as usize
    }
    fn new(val: usize) -> Self {
        AtomicU64::new(val as u64)
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

impl SingleColoredKmers {

    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    // This is used to identify files serialized from this struct
    fn serialization_magic_constant() -> [u8; 4] {
        [54, 229, 250, 84] // Four randomly generated bytes
    }

    // This is used to identify the version of the serialization format
    fn serialization_version_number() -> u32 {
        1_u32
    }

    pub fn serialize(&self, mut out: &mut impl Write) {
        out.write_all(&Self::serialization_magic_constant()).unwrap();
        out.write_all(&Self::serialization_version_number().to_le_bytes()).unwrap();

        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();

        bincode::serialize_into(&mut out, &self.colors).unwrap();
        bincode::serialize_into(&mut out, &self.n_colors).unwrap();
    }

    pub fn load(mut input: &mut impl Read) -> SingleColoredKmers {

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
        let lcs = sbwt::LcsArray::load(input).unwrap();

        let colors = bincode::deserialize_from(&mut input).unwrap();
        let n_colors = bincode::deserialize_from(&mut input).unwrap();

        SingleColoredKmers{sbwt, lcs, colors, n_colors}
    }

    #[allow(dead_code)]
    pub fn n_colors(&self) -> usize {
        self.n_colors
    }
    
    pub fn get_color(&self, colex: usize) -> ColorVecValue {
        assert!(colex < self.sbwt.n_sets());
        self.colors.get_color(colex)
    }

    // Returns an iterator giving the color of each of the n-k+1 k-mers of the query.
    // If query is shorter than k, returns an empty.
    pub fn lookup_kmers<'a, 'b>(&'a self, query: &'b [u8]) -> KmerLookupIterator<'a, 'b>{
        let k = self.sbwt.k();
        let si = StreamingIndex::new(&self.sbwt, &self.lcs);

        let mut ms_iter = si.matching_statistics_iter(query);

        // Skip over the first k-1 positions
        for _ in 0..k-1 {
            ms_iter.next(); // If the iterator ends early, will keep returning None
        }
        KmerLookupIterator { matching_stats_iter: ms_iter, index: self }

    }

    fn progress_print_thread(n_bases_processed: &AtomicU64, quit_signal: Receiver<bool>) {
        log::info!("Processing up to 2n bases, where n is the number of bases in the input."); // 2n due to reverse complements
        let print_interval = 10; // seconds
        let start_time = Instant::now();
        loop {
            match quit_signal.recv_timeout(Duration::from_secs(print_interval)) {
                Ok(_) => return, // Received quit signal
                Err(RecvTimeoutError::Timeout) => { // print_interval seconds has passed
                    let count = n_bases_processed.load(std::sync::atomic::Ordering::Relaxed);
                    let elapsed = start_time.elapsed().as_secs_f64();
                    log::info!("{} bases processed ({:.2} Mbases/sec)", count, count as f64 / elapsed / 1e6);
                },
                Err(RecvTimeoutError::Disconnected) => {
                    // I'm not sure when this would happen, but let's just quit
                    return 
                }
            }
        }
    }

    // Generic function that works on any of u8, 16, u32 and u64
    fn mark_colors<T: SeqStream + Send, A: AtomicColorVec + Send + Sync>(sbwt: &sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: &sbwt::LcsArray,  input_streams: Vec<T>, n_threads: usize) -> ColorStorage {

        let color_ids = A::new(sbwt.n_sets());
        let si = StreamingIndex::new(sbwt, lcs);
        let k = sbwt.k();
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
                    SingleColoredKmers::progress_print_thread(n_bases_processed, quit_print_recv);
                }
            });

            input_streams.into_iter().enumerate().par_bridge().map(|(color, mut input_stream)| {
                let mut thread_progress = 0_usize;
                let mut store_count = 0_usize;
                while let Some(seq) = input_stream.stream_next() {
                    // Figure out colex ranks of k-mers in this sequence
                    let ms = si.matching_statistics_iter(seq);
                    ms.enumerate().for_each(|(i, (len, range))| { 
                        if len == k {
                            debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
                            color_ids.update(range.start, color);
                            store_count += 1;
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
                            n_bases_processed.fetch_add(10000, std::sync::atomic::Ordering::Relaxed);
                            thread_progress = 0;
                        }
                    });
                }

                // Record remaining progress count
                n_bases_processed.fetch_add(thread_progress as u64, std::sync::atomic::Ordering::Relaxed);
                store_count 
            }).sum::<usize>();

            quit_print_send.send(true).unwrap(); // Tell the progress printer to quit (otherwise we hang)
        })});

        // Compress color_ids into a BitVec
        log::info!("Bitpacking color id array");
        let mut compressed_colors = ColorStorage::new(sbwt.n_sets(), n_colors);
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

        let required_bit_width = ColorStorage::required_bit_width(n_colors);

        log::info!("Marking colors");
        let color_storage = if required_bit_width <= 8 {
            SingleColoredKmers::mark_colors::<T, Vec::<AtomicU8>>(&sbwt, &lcs, input_streams, n_threads)
        } else if required_bit_width <= 16 {
            SingleColoredKmers::mark_colors::<T, Vec::<AtomicU16>>(&sbwt, &lcs, input_streams, n_threads)
        } else if required_bit_width <= 32 {
            SingleColoredKmers::mark_colors::<T, Vec::<AtomicU32>>(&sbwt, &lcs, input_streams, n_threads)
        } else {
            SingleColoredKmers::mark_colors::<T, Vec::<AtomicU64>>(&sbwt, &lcs, input_streams, n_threads)
        };

        SingleColoredKmers{
            sbwt, lcs, n_colors, colors: color_storage
        }
    }
}

fn log2_ceil(x: usize) -> usize {
    assert!(x > 0);
    let mut v = 1_usize;
    let mut bits = 0_usize;
    while v < x {
        v <<= 1;
        bits += 1;
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log2_ceil(){
        assert_eq!(log2_ceil(1), 0);
        assert_eq!(log2_ceil(2), 1);
        assert_eq!(log2_ceil(3), 2);
        assert_eq!(log2_ceil(4), 2);
        assert_eq!(log2_ceil(5), 3);
        assert_eq!(log2_ceil(6), 3);
        assert_eq!(log2_ceil(7), 3);
        assert_eq!(log2_ceil(8), 3);

    }
}