use std::ops::Range;

use bitvec_sds::{rank_support_v::RankSupportV, traits::{Pat1, RankSupport, SelectSupport}, wavelet_tree::SelectSupportBoth};

/// Wrapper for a bitvec_sds Wavelet tree. We need to wrap it so that we can
/// implement the foreign ContracLeft trait for it.
#[derive(Debug, Clone)]
pub struct WaveletTree {
    inner: bitvec_sds::wavelet_tree::WaveletTree<RankSupportV<Pat1>, SelectSupportBoth>
}

impl WaveletTree {
    pub fn new(elements: &[u32], n_values_supported: usize) -> Self {
        let inner = bitvec_sds::wavelet_tree::WaveletTree::<RankSupportV::<Pat1>, SelectSupportBoth>::new(&elements, 0, n_values_supported as u32,
            RankSupportV::new,
            SelectSupportBoth::new
        );
        Self { inner }
    }

    pub fn serialize(&self, mut writer: &mut impl std::io::Write) {
        self.inner.serialize(&mut writer);
    }

    pub fn load(mut reader: &mut impl std::io::Read) -> Self {
        let inner = bitvec_sds::wavelet_tree::WaveletTree::<RankSupportV::<Pat1>, SelectSupportBoth>::load(&mut reader);
        Self { inner }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn next_smaller(&self, i: usize, x: u32) -> Option<usize> {
        self.inner.next_smaller(i, x)
    }

    pub fn prev_smaller(&self, i: usize, x: u32) -> Option<usize> {
        self.inner.prev_smaller(i, x)
    }

    pub fn access(&self, i: usize) -> u32 {
        self.inner.access(i)
    }

    pub fn range_min(&self, l: usize, r: usize) -> Option<u32> {
        self.inner.range_min(l, r)
    }

    pub fn rank(&self, x: u32, i: usize) -> usize {
        self.inner.rank(x, i)
    }
    pub fn range_rank(&self, x: u32, l: usize, r: usize) -> usize {
        self.inner.range_rank(x, l, r)
    }

    pub fn value_range(&self) -> Range<usize> {
        self.inner.value_range()
    }
}
impl sbwt::ContractLeft for WaveletTree {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {

        let new_start = match self.inner.prev_smaller(I.start + 1, target_len as u32) {
            None => 0,
            Some(s) => s,
        };

        assert!(I.end > 0);
        let new_end = match self.inner.next_smaller(I.end - 1, target_len as u32) {
            None => self.inner.len(),
            Some(s) => s,
        };

        new_start..new_end
    }
}


