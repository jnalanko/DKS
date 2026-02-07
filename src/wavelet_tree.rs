
/*
impl<R: RankSupport, S: SelectSupport> sbwt::ContractLeft for WaveletTree<R,S> {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {

        let new_start = match self.prev_smaller(I.start + 1, target_len as u32) {
            None => 0,
            Some(s) => s,
        };

        assert!(I.end > 0);
        let new_end = match self.next_smaller(I.end - 1, target_len as u32) {
            None => self.len(),
            Some(s) => s,
        };

        new_start..new_end
    }
}
*/


