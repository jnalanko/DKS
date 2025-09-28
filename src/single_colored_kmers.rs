use std::io::{Read, Write};

use bitvec::prelude::*;
use bitvec::{field::BitField, order::Lsb0, vec::BitVec};
use sbwt::{SbwtIndex, SeqStream, StreamingIndex};
use serde::Serialize;

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

impl SingleColoredKmers {

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
    
    pub fn get_color(&self, colex: usize) -> usize {
        assert!(colex < self.sbwt.n_sets());
        self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].load_le()
    }

    pub fn new<T: SeqStream>(sbwt: sbwt::SbwtIndex<sbwt::SubsetMatrix>, lcs: sbwt::LcsArray, input_streams: Vec<T>) -> Self {
        let si = StreamingIndex::new(&sbwt, &lcs);
        let k = sbwt.k();
        let n_colors = input_streams.len();
        let bits_per_color = (n_colors.next_power_of_two().trailing_zeros()+1) as usize;
        let mut color_ids: BitVec::<usize, Lsb0> = bitvec![0; sbwt.n_sets()*bits_per_color];
        for (color, mut input_stream) in input_streams.into_iter().enumerate() {
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
                    color_ids[colex*bits_per_color .. (colex+1)*bits_per_color].store_le(color);
                }
            }
        }

        SingleColoredKmers{
            sbwt, lcs, n_colors, bits_per_color, colors: color_ids
        }
    }
}