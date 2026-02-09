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
    n_colors: usize, // Number of single colors, not including the special colors for multiple and none.
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
    
    fn get_color_of_range(&self, range: Range<usize>) -> ColorVecValue {
        // This is O(|range| in the worst case)

        if range.is_empty() { return ColorVecValue::None }

        let mut existing_color: Option<usize> = None;
        for colex in range {
            match self.get_color(colex) {
                ColorVecValue::Single(c) => {
                    match existing_color {
                        Some(e) => {
                            if c != e {
                                return ColorVecValue::Multiple;
                            }    
                        },
                        None => existing_color = Some(c),
                    }
                },
                ColorVecValue::Multiple => return ColorVecValue::Multiple,
                ColorVecValue::None => (),
            }
        }

        match existing_color {
            Some(c) => ColorVecValue::Single(c),
            None => ColorVecValue::None,
        }
    }
}

impl SimpleColorStorage {

    pub fn set_color(&mut self, colex: usize, color: ColorVecValue) {
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

    pub fn new(len: usize, n_colors: usize) -> Self {
        let bits_per_color = Self::required_bit_width(n_colors);
        SimpleColorStorage {
            n_colors,
            colors: bitvec![0; len * bits_per_color],
            bits_per_color,
        }
    }

    pub fn required_bit_width(n_colors: usize) -> usize {
        log2_ceil(n_colors + 2) // +2 to reserve two special values: one for "no color" and one for "multiple colors"
    }

    pub fn into_wavelet_tree_storage(self) -> WTColorStorage {
        assert!(self.bits_per_color <= 32);
        let value_range_end = 1 << self.bits_per_color; // Exclusive end of the supported value range

        let wt = crate::wavelet_tree::WaveletTreeWrapper::new(self, value_range_end);
        WTColorStorage { colors: wt }
    }

    pub fn n_colors(&self) -> usize {
        self.n_colors
    }
}

#[derive(Debug, Clone)]
pub struct WTColorStorage {
    colors: crate::wavelet_tree::WaveletTreeWrapper,
}

impl WTColorStorage {
    fn get_color(&self, colex: usize) -> ColorVecValue {
        let x = self.colors.access(colex) as usize;
        let none_id = self.colors.value_range().end-1; // Max value is reserved for None
        let multiple_id = none_id - 1; // Max - 1 is reserved for Multiple
        if x == none_id { 
            ColorVecValue::None
        } else if x == multiple_id { 
            ColorVecValue::Multiple
        } else {
            ColorVecValue::Single(x)
        }
    }

    fn get_color_of_range(&self, range: Range<usize>) -> ColorVecValue {
        if range.is_empty() { return ColorVecValue::None }

        let none_id = self.colors.value_range().end-1; // Max value is reserved for None
        let multiple_id = none_id - 1; // Max - 1 is reserved for Multiple

        let (a,b) = self.colors.range_bottom2(range.start, range.end);

        match (a,b) {
            (None, None) => {
                assert!(range.is_empty());
                ColorVecValue::None
            }
            (None, Some(_)) => panic!("range_bottom2 returned (None, Some)"), // Should never happen
            (Some(x), None) => {
                let x = x as usize;
                if x == none_id { ColorVecValue::None }
                else if x == multiple_id { ColorVecValue::Multiple }
                else { ColorVecValue::Single(x) }
            },
            (Some(x), Some(y)) => {
                let x = x as usize;
                let y = y as usize; // If x or y is none_id, it's this, because none_id is the largest possible
                if y == none_id && x != multiple_id {
                    ColorVecValue::Single(x)
                } else {
                    ColorVecValue::Multiple
                }
            },
        }
    }

}

impl ColorStorage for WTColorStorage {
    fn get_color(&self, colex: usize) -> ColorVecValue {
        self.get_color(colex)
    }

    fn get_color_of_range(&self, range: Range<usize>) -> ColorVecValue {
        self.get_color_of_range(range)
    }
}

impl MySerialize for WTColorStorage {
    fn serialize(&self, out: &mut impl Write) {
        self.colors.serialize(out);
    }

    fn load(input: &mut impl Read) -> Box<Self> {
        let colors = WaveletTreeWrapper::load(input);
        Box::new(WTColorStorage { colors })
    }
}

impl From<SimpleColorStorage> for WTColorStorage {
    fn from(s: SimpleColorStorage) -> Self {
        s.into_wavelet_tree_storage() 
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