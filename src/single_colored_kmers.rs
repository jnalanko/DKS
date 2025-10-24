use std::io::{Read, Write};
use std::sync::{Arc, Mutex};

use bitvec::prelude::*;
use bitvec::{field::BitField, order::Lsb0, vec::BitVec};
use sbwt::{LcsArray, MatchingStatisticsIterator, SbwtIndex, SeqStream, StreamingIndex, SubsetMatrix};
use rayon::prelude::*;

// This bit vector of length 256 marks the ascii values of these characters: acgtACGT
const IS_DNA: BitArray<[u32; 8]> = bitarr![const u32, Lsb0; 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

#[derive(Debug, Clone)]
pub struct SingleColoredKmers {
    sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>,
    lcs: sbwt::LcsArray,
    colors: BitVec::<usize, Lsb0>,
    n_colors: usize,
    bits_per_color: usize,
}

pub struct KmerLookupIterator<'a, 'b> {
    // This iterator should be initialized so that the first k-1 MS values are skipped
    matching_stats_iter: MatchingStatisticsIterator<'a, 'b, SbwtIndex::<SubsetMatrix>, LcsArray>,
    index: &'a SingleColoredKmers,
}

impl Iterator for KmerLookupIterator<'_, '_> {
    type Item = Option<usize>; // Color id of k-mer, if exists

    fn next(&mut self) -> Option<Self::Item> {
        let (len, range) = self.matching_stats_iter.next()?;

        if len == self.index.sbwt.k() {
            debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
            let colex = range.start;
            Some(Some(self.index.get_color(colex)))
        } else {
            Some(None) // Iterator not finished but the k-mer is not found
        }
    }
}

impl SingleColoredKmers {

    pub fn k(&self) -> usize {
        self.sbwt.k()
    }

    pub fn serialize(&self, mut out: &mut impl Write) {
        self.sbwt.serialize(out).unwrap();
        self.lcs.serialize(out).unwrap();

        // Todo: specify format (is now: usize, should be: little endian 64 bit)
        // Todo: magic string to check the file format.
        bincode::serialize_into(&mut out, &self.colors).unwrap();
        bincode::serialize_into(&mut out, &self.n_colors).unwrap();
        bincode::serialize_into(&mut out, &self.bits_per_color).unwrap();
    }

    pub fn load(mut input: &mut impl Read) -> SingleColoredKmers {
        let sbwt = SbwtIndex::<sbwt::SubsetMatrix>::load(input).unwrap();
        let lcs = sbwt::LcsArray::load(input).unwrap();

        let colors = bincode::deserialize_from(&mut input).unwrap();
        let n_colors = bincode::deserialize_from(&mut input).unwrap();
        let bits_per_color = bincode::deserialize_from(&mut input).unwrap();

        SingleColoredKmers{sbwt, lcs, colors, n_colors, bits_per_color}
    }

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }
    
    pub fn get_color(&self, colex: usize) -> usize {
        assert!(colex < self.sbwt.n_sets());
        self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].load_le()
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

    fn write_ints(bv: &mut BitVec::<usize, Lsb0>, indices: &[usize], value: usize, bit_width: usize) {
        for &i in indices.iter() {
            bv[i*bit_width..(i+1)*bit_width].store_le(value);
        }
    }

    pub fn new<T: SeqStream + Send>(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, input_streams: Vec<T>, n_threads: usize) -> Self {

        let si = StreamingIndex::new(&sbwt, &lcs);
        let k = sbwt.k();
        let n_colors = input_streams.len();
        let bits_per_color = (n_colors.next_power_of_two().trailing_zeros()+1) as usize;
        let color_ids: BitVec::<usize, Lsb0> = bitvec![0; sbwt.n_sets()*bits_per_color];
        let color_ids = Arc::new(Mutex::new(color_ids)); // Wrap in a mutex for access from different threads

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
        thread_pool.install(|| {
            input_streams.into_iter().enumerate().par_bridge().for_each(|(color,mut input_stream)| {
                let mut store_buffer = Vec::<usize>::new(); // Buffer of writes to the shared color_ids vector
                let store_buffer_cap = 10000; // Flush buffer when capacity reaches this
                while let Some(seq) = input_stream.stream_next() {
                    // Figure out colex ranks of k-mers in this sequence
                    let ms = si.matching_statistics(seq);
                    let colex_iter = ms.iter().enumerate().filter_map(|(i, (len, range))| if 
                        *len == k {
                            debug_assert!(range.len() == 1); // Full k-mer should have a singleton range
                            Some(range.start)
                        } else {
                            if i >= k-1 {
                                // All valid k-mers should be found. If we're here, the k-mer must have had non-ACGT
                                // characters which make it invalid. Let's verify that.
                                let kmer = &seq[i-(k-1)..=i];
                                let all_ACGT = kmer.iter().all(|c| IS_DNA[*c as usize]);
                                if all_ACGT {
                                    panic!("Error: k-mer {} not found in sbwt", String::from_utf8_lossy(&kmer));
                                }
                            }
                            None
                        }
                    );

                    // Store color ids
                    for colex in colex_iter {
                        store_buffer.push(colex);
                        if store_buffer.len() >= store_buffer_cap {
                            // Flush buffer
                            Self::write_ints(&mut color_ids.lock().unwrap(), &store_buffer, color, bits_per_color);
                            store_buffer.clear();
                        }
                    }
                }

                // Flush the remaining writes in the buffer
                Self::write_ints(&mut color_ids.lock().unwrap(), &store_buffer, color, bits_per_color);
            })
        });

        // Only one reference should exists now because all the threads are finished
        let color_ids = Arc::into_inner(color_ids).unwrap().into_inner().unwrap(); 
        SingleColoredKmers{
            sbwt, lcs, n_colors, bits_per_color, colors: color_ids
        }
    }
}