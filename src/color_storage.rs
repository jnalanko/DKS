use crate::lca_tree::LcaTree;
use crate::traits::*;
use crate::wavelet_tree::WaveletTreeWrapper;
use serde::{Serialize, Deserialize};
use bitvec::prelude::*;
use std::io::{Read, Write};
use std::ops::Range;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleColorStorage {
    colors: BitVec::<usize, Lsb0>,
    bits_per_color: usize,
    n_colors: usize, // Number of colors, not including the special "none" color.
}

impl bitvec_sds::traits::RandomAccessU32 for SimpleColorStorage {
    fn len(&self) -> usize {
        self.colors.len() / self.bits_per_color
    }

    fn get(&self, idx: usize) -> u32 {
        self.colors[idx*self.bits_per_color .. (idx+1)*self.bits_per_color].load_le()
    }
}

impl MySerialize for SimpleColorStorage {
    fn serialize(&self, mut out: &mut impl Write) {
        bincode::serialize_into(&mut out, &self.n_colors).unwrap();
        bincode::serialize_into(&mut out, &self.bits_per_color).unwrap();
        bincode::serialize_into(&mut out, &self.colors).unwrap();
    }

    fn load(mut input: &mut impl Read) -> Box<Self> {
        let n_colors: usize = bincode::deserialize_from(&mut input).unwrap();
        let bits_per_color: usize = bincode::deserialize_from(&mut input).unwrap();
        let colors: BitVec<usize, Lsb0> = bincode::deserialize_from(&mut input).unwrap();
        Box::new(SimpleColorStorage { n_colors, colors, bits_per_color })
    }
}

impl ColorStorage for SimpleColorStorage {
    fn get_color(&self, colex: usize) -> Option<usize> {
        let x: usize = self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].load_le();
        if x == (1 << self.bits_per_color) - 1 { // Max value is reserved for None
            None
        } else {
            Some(x)
        }
    }

    fn set_color(&mut self, colex: usize, value: Option<usize>) {
        let x = match value {
            None => (1 << self.bits_per_color) - 1,
            Some(x) => {
                assert!(x < (1 << self.bits_per_color) - 1);
                x
            }
        };
        
        self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].store_le(x);
    }
    
    fn get_color_of_range(&self, range: Range<usize>, color_hierarchy: &LcaTree) -> Option<usize> {
        // This is O(|range| in the worst case)

        if range.is_empty() { return None }

        let mut lca: Option<usize> = None;
        for colex in range {
            lca = color_hierarchy.lca_options(lca, self.get_color(colex));
            if lca == Some(color_hierarchy.root()) {
                // Already at root -> can early exit here
                return lca
            }
        }
        lca
    }
}

impl SimpleColorStorage {

    pub fn set_color(&mut self, colex: usize, color: Option<usize>) {
        let value = match color {
            None => (1 << self.bits_per_color) - 1, // Max value is reserved for None
            Some(x) => {
                assert!(x < (1 << self.bits_per_color) - 1);
                x
            }
        };
        self.colors[colex*self.bits_per_color .. (colex+1)*self.bits_per_color].store_le(value);
    }

    pub fn new(len: usize, n_colors: usize) -> Self {
        let bits_per_color = Self::required_bit_width(n_colors);
        SimpleColorStorage {
            n_colors,
            colors: bitvec![0; len * bits_per_color],
            bits_per_color,
        }
    }

    pub fn required_bit_width(n_colors: usize) -> usize {
        log2_ceil(n_colors + 1) // +1 is for the special "none" value 
    }

    pub fn n_colors(&self) -> usize {
        self.n_colors
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