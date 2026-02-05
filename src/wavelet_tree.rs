use std::sync::Arc;

use bitvec::prelude::*;
use bps_sada::traits::*;
use bps_sada::rank_support_v::RankSupportV;
use bps_sada::select_support_mcl::{self, SelectSupportMcl};

/// Rank support: # of 1-bits in B[0..idx)
pub trait RankSupport {
    fn rank1(&self, idx: usize) -> usize;

    #[inline]
    fn rank0(&self, idx: usize) -> usize {
        idx - self.rank1(idx)
    }
}

impl RankSupport for RankSupportV<Pat1> {
    fn rank1(&self, idx: usize) -> usize {
        self.rank(idx)
    }
}

/// Select support: position of k-th bit (1-based k), returning 0-based position.
pub trait SelectSupport {
    fn select1(&self, k: usize) -> Option<usize>;
    fn select0(&self, k: usize) -> Option<usize>;
}

struct SelectSupportMcl0and1 {
    sel0: SelectSupportMcl<Sel0>,
    sel1: SelectSupportMcl<Sel1>,
}

impl SelectSupport for SelectSupportMcl0and1 {
    fn select1(&self, k: usize) -> Option<usize> {
        Some(self.sel1.select(k)) // Will panic if k is larger than the number of ones
    }

    fn select0(&self, k: usize) -> Option<usize> {
        Some(self.sel0.select(k)) // Will panic if k is larger than the number of zeros
    }
}


struct Node<R, S> {
    lo: u32,
    hi: u32, // exclusive
    mid: u32,
    bits: Arc<BitVec<u64, Lsb0>>,
    rank: R,
    sel: S,
    left: Option<usize>,
    right: Option<usize>,
}

pub struct WaveletTree<R, S> {
    nodes: Vec<Node<R, S>>,
    n: usize,
    lo: u32,
    hi: u32, // exclusive
}

impl<R, S> WaveletTree<R, S>
where
    R: RankSupport,
    S: SelectSupport,
{
    /// Build a wavelet tree for values in [lo, hi).
    ///
    /// `build_rank` and `build_sel` construct rank/select structures for each node's bitvector.
    pub fn new<FR, FS>(data: &[u32], lo: u32, hi: u32, mut build_rank: FR, mut build_sel: FS) -> Self
    where
        FR: FnMut(Arc<BitVec<u64, Lsb0>>) -> R,
        FS: FnMut(Arc<BitVec<u64, Lsb0>>) -> S,
    {
        assert!(lo < hi, "alphabet range must be non-empty");
        for &v in data {
            assert!(v >= lo && v < hi, "value {v} out of range [{lo},{hi})");
        }

        let mut wt = WaveletTree {
            nodes: Vec::new(),
            n: data.len(),
            lo,
            hi,
        };

        fn build_rec<R, S, FR, FS>(
            wt: &mut WaveletTree<R, S>,
            seq: &[u32],
            lo: u32,
            hi: u32,
            build_rank: &mut FR,
            build_sel: &mut FS,
        ) -> usize
        where
            R: RankSupport,
            S: SelectSupport,
            FR: FnMut(Arc<BitVec<u64, Lsb0>>) -> R,
            FS: FnMut(Arc<BitVec<u64, Lsb0>>) -> S,
        {
            let mid = lo + (hi - lo) / 2;

            // Leaf interval size 1.
            if hi - lo == 1 {
                let bits = Arc::new(BitVec::<u64, Lsb0>::new());
                let rank = build_rank(bits.clone());
                let sel = build_sel(bits.clone());
                wt.nodes.push(Node {
                    lo,
                    hi,
                    mid,
                    bits,
                    rank,
                    sel,
                    left: None,
                    right: None,
                });
                return wt.nodes.len() - 1;
            }

            let mut bits = BitVec::<u64, Lsb0>::with_capacity(seq.len());
            let mut left_vals = Vec::new();
            let mut right_vals = Vec::new();
            left_vals.reserve(seq.len());
            right_vals.reserve(seq.len());

            for &v in seq {
                let go_right = v >= mid;
                bits.push(go_right);
                if go_right {
                    right_vals.push(v);
                } else {
                    left_vals.push(v);
                }
            }

            let bits = Arc::new(bits);

            let rank = build_rank(bits.clone());
            let sel = build_sel(bits.clone());

            let node_idx = wt.nodes.len();
            wt.nodes.push(Node {
                lo,
                hi,
                mid,
                bits,
                rank,
                sel,
                left: None,
                right: None,
            });

            let left = build_rec(wt, &left_vals, lo, mid, build_rank, build_sel);
            let right = build_rec(wt, &right_vals, mid, hi, build_rank, build_sel);
            wt.nodes[node_idx].left = Some(left);
            wt.nodes[node_idx].right = Some(right);

            node_idx
        }

        let _root = build_rec(&mut wt, data, lo, hi, &mut build_rank, &mut build_sel);
        wt
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    /// Next position > i with value < x in O(log k).
    pub fn next_smaller(&self, i: usize, x: u32) -> Option<usize> {
        if self.n == 0 || i + 1 >= self.n {
            return None;
        }
        if x <= self.lo {
            return None;
        }
        // Search in suffix [i+1, n)
        self.next_in_value_prefix_range(0, i + 1, self.n, x)
    }

    /// Previous position < i with value < x in O(log k).
    pub fn prev_smaller(&self, i: usize, x: u32) -> Option<usize> {
        if self.n == 0 || i == 0 {
            return None;
        }
        if x <= self.lo {
            return None;
        }
        // Search in prefix [0, i)
        self.prev_in_value_prefix_range(0, 0, i, x)
    }

    fn next_in_value_prefix_range(
        &self,
        node_idx: usize,
        l: usize,
        r: usize,
        x: u32,
    ) -> Option<usize> {
        if l >= r {
            return None;
        }

        let node = &self.nodes[node_idx];

        if x <= node.lo {
            return None;
        }
        if x >= node.hi {
            // FULL value coverage => earliest position in this node interval
            return Some(l);
        }

        if node.hi - node.lo == 1 {
            // x > node.lo already holds, so anything qualifies; earliest is l
            return Some(l);
        }

        let left_idx = node.left.expect("internal node must have left");
        let right_idx = node.right.expect("internal node must have right");

        let l0 = node.rank.rank0(l);
        let r0 = node.rank.rank0(r);
        let l1 = node.rank.rank1(l);
        let r1 = node.rank.rank1(r);

        let cand_left = self
            .next_in_value_prefix_range(left_idx, l0, r0, x)
            .and_then(|p_child| node.sel.select0(p_child + 1));

        let cand_right = self
            .next_in_value_prefix_range(right_idx, l1, r1, x)
            .and_then(|p_child| node.sel.select1(p_child + 1));

        match (cand_left, cand_right) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (Some(a), None) => Some(a),
            (None, Some(b)) => Some(b),
            (None, None) => None,
        }
    }

    fn prev_in_value_prefix_range(
        &self,
        node_idx: usize,
        l: usize,
        r: usize,
        x: u32,
        ) -> Option<usize> {
            if l >= r {
                return None;
            }

            let node = &self.nodes[node_idx];

            if x <= node.lo {
                return None;
            }
            if x >= node.hi {
                // FULL value coverage => latest position in this node interval
                return Some(r - 1);
            }

            if node.hi - node.lo == 1 {
                // x > node.lo holds => anything qualifies; latest is r-1
                return Some(r - 1);
            }

            let left_idx = node.left.expect("internal node must have left");
            let right_idx = node.right.expect("internal node must have right");

            let l0 = node.rank.rank0(l);
            let r0 = node.rank.rank0(r);
            let l1 = node.rank.rank1(l);
            let r1 = node.rank.rank1(r);

            let cand_left = self
                .prev_in_value_prefix_range(left_idx, l0, r0, x)
                .and_then(|p_child| node.sel.select0(p_child + 1));

            let cand_right = self
                .prev_in_value_prefix_range(right_idx, l1, r1, x)
                .and_then(|p_child| node.sel.select1(p_child + 1));

            match (cand_left, cand_right) {
                (Some(a), Some(b)) => Some(a.max(b)),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            }
        }

    /// Return the first position in [l, r) within this node (mapped to this node's coords).
    /// O(log k) via select mapping down/up.
    fn first_in_interval(&self, node_idx: usize, l: usize, r: usize) -> usize {
        let node = &self.nodes[node_idx];
        if node.hi - node.lo == 1 {
            return l;
        }

        let left_idx = node.left.expect("internal node must have left");
        let right_idx = node.right.expect("internal node must have right");

        let l0 = node.rank.rank0(l);
        let r0 = node.rank.rank0(r);
        if l0 < r0 {
            let p_child = self.first_in_interval(left_idx, l0, r0);
            node.sel
                .select0(p_child + 1)
                .expect("select0 must exist for valid child position")
        } else {
            let l1 = node.rank.rank1(l);
            let r1 = node.rank.rank1(r);
            // Since [l,r) non-empty and left empty, right must be non-empty.
            let p_child = self.first_in_interval(right_idx, l1, r1);
            node.sel
                .select1(p_child + 1)
                .expect("select1 must exist for valid child position")
        }
    }

    /// Return the last position in [l, r) within this node (mapped to this node's coords).
    fn last_in_interval(&self, node_idx: usize, l: usize, r: usize) -> usize {
        let node = &self.nodes[node_idx];
        if node.hi - node.lo == 1 {
            return r - 1;
        }

        let left_idx = node.left.expect("internal node must have left");
        let right_idx = node.right.expect("internal node must have right");

        let l1 = node.rank.rank1(l);
        let r1 = node.rank.rank1(r);
        if l1 < r1 {
            let p_child = self.last_in_interval(right_idx, l1, r1);
            node.sel
                .select1(p_child + 1)
                .expect("select1 must exist for valid child position")
        } else {
            let l0 = node.rank.rank0(l);
            let r0 = node.rank.rank0(r);
            let p_child = self.last_in_interval(left_idx, l0, r0);
            node.sel
                .select0(p_child + 1)
                .expect("select0 must exist for valid child position")
        }
    }
}

#[cfg(test)]
mod stress {
    use std::sync::Arc;

    use super::*;

    // --- Deterministic RNG (fast, no deps) ---

    #[derive(Clone)]
    struct SplitMix64(u64);
    impl SplitMix64 {
        fn new(seed: u64) -> Self {
            Self(seed)
        }
        fn next_u64(&mut self) -> u64 {
            // splitmix64
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
        fn next_usize(&mut self) -> usize {
            (self.next_u64() >> 1) as usize
        }
        fn gen_range_u32(&mut self, lo: u32, hi: u32) -> u32 {
            assert!(lo < hi);
            let span = (hi - lo) as u64;
            lo + (self.next_u64() % span) as u32
        }
        fn gen_range_usize(&mut self, lo: usize, hi: usize) -> usize {
            assert!(lo < hi);
            let span = (hi - lo) as u64;
            lo + (self.next_u64() % span) as usize
        }
        fn coin(&mut self, num: u64, den: u64) -> bool {
            debug_assert!(num <= den && den > 0);
            (self.next_u64() % den) < num
        }
    }

    // --- Brute references ---

    fn brute_next_smaller(a: &[u32], i: usize, x: u32) -> Option<usize> {
        if i + 1 >= a.len() {
            return None;
        }
        for j in (i + 1)..a.len() {
            if a[j] < x {
                return Some(j);
            }
        }
        None
    }

    fn brute_prev_smaller(a: &[u32], i: usize, x: u32) -> Option<usize> {
        if i == 0 {
            return None;
        }
        for j in (0..i).rev() {
            if a[j] < x {
                return Some(j);
            }
        }
        None
    }

    // --- Data generators that are intentionally nasty for wavelet trees ---

    /// Generate data in [0, k) with a mixture of:
    /// - long runs (highly compressible)
    /// - alternating patterns (maximizes crossings between partitions)
    /// - heavy skew (most values near 0 or k-1)
    /// - occasional uniform noise
    fn gen_adversarial_array(rng: &mut SplitMix64, n: usize, k: u32) -> Vec<u32> {
        assert!(k >= 2);
        let mut a = Vec::with_capacity(n);
        let mut mode = 0u8;

        let mut i = 0usize;
        while i < n {
            // Switch mode sometimes
            if rng.coin(1, 20) {
                mode = (rng.next_u64() % 5) as u8;
            }

            match mode {
                // 0: long run of same value
                0 => {
                    let v = if rng.coin(1, 2) { 0 } else { k - 1 };
                    let run = rng.gen_range_usize(1, (n - i).min(512) + 1);
                    for _ in 0..run {
                        a.push(v);
                    }
                    i += run;
                }
                // 1: alternating low/high
                1 => {
                    let run = rng.gen_range_usize(1, (n - i).min(1024) + 1);
                    for t in 0..run {
                        a.push(if (t & 1) == 0 { 0 } else { k - 1 });
                    }
                    i += run;
                }
                // 2: staircase: 0,0,1,1,2,2,... wrapping
                2 => {
                    let run = rng.gen_range_usize(1, (n - i).min(1024) + 1);
                    let mut v = rng.gen_range_u32(0, k);
                    for t in 0..run {
                        // change every 2 steps
                        if (t % 2) == 0 {
                            v = (v + 1) % k;
                        }
                        a.push(v);
                    }
                    i += run;
                }
                // 3: heavy skew to small numbers, occasional spikes
                3 => {
                    let run = rng.gen_range_usize(1, (n - i).min(1024) + 1);
                    for _ in 0..run {
                        let r = rng.next_u64() % 100;
                        let v = if r < 85 {
                            // very small
                            (rng.next_u64() % 8) as u32 % k
                        } else if r < 95 {
                            // very large
                            k - 1 - ((rng.next_u64() % 8) as u32 % k)
                        } else {
                            // uniform
                            rng.gen_range_u32(0, k)
                        };
                        a.push(v);
                    }
                    i += run;
                }
                // 4: uniform noise
                _ => {
                    let run = rng.gen_range_usize(1, (n - i).min(2048) + 1);
                    for _ in 0..run {
                        a.push(rng.gen_range_u32(0, k));
                    }
                    i += run;
                }
            }
        }

        a
    }

    // --- Stress test ---

    /// A much more challenging stress test:
    /// - Larger n
    /// - Larger k (log k bigger)
    /// - Adversarial distributions
    /// - Many random queries, including edge and boundary x values
    /// - Some targeted x values around midpoints and around present values
    #[test]
    fn stress_next_prev_smaller_hard() {
        let mut rng = SplitMix64::new(0xD1CE_BA5E_F00D_F00Du64);

        // Keep this reasonably sized for unit tests; bump n/iters locally if desired.
        // If you want a "burn-in" stress test, set n=200_000, queries=200_000.
        let n: usize = 50_000;
        let k: u32 = 1 << 16; // 65536 alphabet -> log2(k)=16 levels

        let a = gen_adversarial_array(&mut rng, n, k);

        let wt = WaveletTree::new(
            &a, 0, k, 
            RankSupportV::<Pat1>::new, 
            |bv| SelectSupportMcl0and1 {
                sel0: SelectSupportMcl::new(bv.clone()),
                sel1: SelectSupportMcl::new(bv.clone())
            }
        );

        // Query budget: mix of random and targeted
        let queries = 80_000usize;

        for t in 0..queries {
            // Choose i with bias toward edges and random interior
            let i = if t % 10 == 0 {
                // edge-heavy
                match t % 4 {
                    0 => 0,
                    1 => n / 2,
                    2 => n - 1,
                    _ => rng.gen_range_usize(0, n),
                }
            } else {
                rng.gen_range_usize(0, n)
            };

            // Choose x: include boundaries, midpoints, near-by values, and random
            let x = match t % 12 {
                0 => 0,
                1 => 1,
                2 => k - 1,
                3 => k,
                4 => k + 1, // outside range on purpose
                5 => k / 2,
                6 => (k / 2).saturating_sub(1),
                7 => (k / 2) + 1,
                8 => a[i].saturating_add(1), // just above current
                9 => a[i],                   // equal to current
                10 => a[i].saturating_sub(1), // just below current
                _ => rng.gen_range_u32(0, k),
            };

            let got_n = wt.next_smaller(i, x);
            let exp_n = brute_next_smaller(&a, i, x);
            assert_eq!(
                got_n, exp_n,
                "next_smaller mismatch @t={}: i={}, x={}, got={:?}, exp={:?}",
                t, i, x, got_n, exp_n
            );

            let got_p = wt.prev_smaller(i, x);
            let exp_p = brute_prev_smaller(&a, i, x);
            assert_eq!(
                got_p, exp_p,
                "prev_smaller mismatch @t={}: i={}, x={}, got={:?}, exp={:?}",
                t, i, x, got_p, exp_p
            );

            // Additional sanity: if an answer exists, verify it is actually < x and in correct direction
            if let Some(j) = got_n {
                assert!(j > i);
                assert!(a[j] < x);
                // ensure it's the next such position
                for jj in (i + 1)..j {
                    assert!(
                        !(a[jj] < x),
                        "next_smaller not minimal: found earlier jj={} with a[jj]={} < x={}",
                        jj,
                        a[jj],
                        x
                    );
                }
            }
            if let Some(j) = got_p {
                assert!(j < i);
                assert!(a[j] < x);
                for jj in (j + 1)..i {
                    assert!(
                        !(a[jj] < x),
                        "prev_smaller not maximal: found later jj={} with a[jj]={} < x={}",
                        jj,
                        a[jj],
                        x
                    );
                }
            }
        }
    }
}

