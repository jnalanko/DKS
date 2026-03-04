use std::io::{self, Read, Write};

/// A rooted tree with nodes labeled 0..n where 0 is the root.
/// Supports LCA queries via a precomputed n^2 lookup table.
/// `parent[0] == 0` (root is its own parent).
#[derive(Debug)]
pub struct LcaTree {
    n: usize,
    /// parent[i] is the parent of node i; parent[0] == 0.
    parent: Vec<usize>,
    /// Flat n×n table: lca_table[i * n + j] = LCA(i, j)
    lca_table: Vec<usize>,
}

impl LcaTree {
    /// Construct a tree from `n` nodes and `edges`, where each edge `(child, parent)`
    /// points toward the root (node 0). Returns an error if the edges do not form a valid tree.
    pub fn new(n: usize, edges: Vec<(usize, usize)>) -> Result<Self, String> {
        if n == 0 {
            return Err("Tree must have at least one node".to_string());
        }

        if edges.len() != n - 1 {
            return Err(format!(
                "Expected {} edges for a tree with {} nodes, got {}",
                n - 1,
                n,
                edges.len()
            ));
        }

        // usize::MAX = "not yet assigned"
        let mut parent = vec![usize::MAX; n];
        parent[0] = 0;
        let mut children = vec![vec![]; n];

        for &(child, par) in &edges {
            if child >= n || par >= n {
                return Err(format!(
                    "Edge ({child}, {par}) contains out-of-range node index (n={n})"
                ));
            }
            if child == par {
                return Err(format!("Self-loop at node {child}"));
            }
            if child == 0 {
                return Err("Root node 0 cannot have a parent".to_string());
            }
            if parent[child] != usize::MAX {
                return Err(format!("Node {child} has multiple parents"));
            }
            parent[child] = par;
            children[par].push(child);
        }

        // BFS from root to compute depths and verify connectivity (no cycles, all reachable).
        let mut depth = vec![usize::MAX; n];
        depth[0] = 0;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(0usize);
        let mut visited = 0usize;

        while let Some(node) = queue.pop_front() {
            visited += 1;
            for &child in &children[node] {
                if depth[child] != usize::MAX {
                    return Err(format!("Cycle detected at node {child}"));
                }
                depth[child] = depth[node] + 1;
                queue.push_back(child);
            }
        }

        if visited != n {
            return Err("Tree is disconnected: not all nodes are reachable from root".to_string());
        }

        // Precompute all LCA pairs.
        let mut lca_table = vec![0usize; n * n];
        for i in 0..n {
            for j in 0..n {
                lca_table[i * n + j] = Self::naive_lca(i, j, &parent, &depth);
            }
        }

        Ok(LcaTree { n, parent, lca_table })
    }

    fn naive_lca(mut a: usize, mut b: usize, parent: &[usize], depth: &[usize]) -> usize {
        while depth[a] > depth[b] {
            a = parent[a];
        }
        while depth[b] > depth[a] {
            b = parent[b];
        }
        while a != b {
            a = parent[a];
            b = parent[b];
        }
        a
    }

    /// Returns the lowest common ancestor of nodes `a` and `b`.
    pub fn lca(&self, a: usize, b: usize) -> usize {
        assert!(a < self.n && b < self.n, "Node index out of bounds");
        self.lca_table[a * self.n + b]
    }

    pub fn n(&self) -> usize {
        self.n
    }

    /// Returns the parent of `node`. For the root, returns 0 (itself).
    pub fn parent(&self, node: usize) -> usize {
        assert!(node < self.n);
        self.parent[node]
    }

    // --- Serialization ---

    /// Write the tree to `w` in a simple binary format:
    /// each `Vec<usize>` is stored as a little-endian u64 length followed by its elements
    /// as little-endian u64 values.
    pub fn serialize<W: Write>(&self, w: &mut W) -> io::Result<()> {
        write_u64(w, self.n as u64)?;
        write_usize_slice(w, &self.parent)?;
        write_usize_slice(w, &self.lca_table)?;
        Ok(())
    }

    /// Load a tree previously written with [`serialize`].
    pub fn load<R: Read>(r: &mut R) -> io::Result<Self> {
        let n = read_u64(r)? as usize;
        let parent = read_usize_vec(r)?;
        let lca_table = read_usize_vec(r)?;
        Ok(LcaTree { n, parent, lca_table })
    }
}

fn write_u64<W: Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn write_usize_slice<W: Write>(w: &mut W, v: &[usize]) -> io::Result<()> {
    write_u64(w, v.len() as u64)?;
    for &x in v {
        write_u64(w, x as u64)?;
    }
    Ok(())
}

fn read_usize_vec<R: Read>(r: &mut R) -> io::Result<Vec<usize>> {
    let len = read_u64(r)? as usize;
    let mut v = Vec::with_capacity(len);
    for _ in 0..len {
        v.push(read_u64(r)? as usize);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: build edges for a path 0-1-2-..-(n-1), edges point toward root (0).
    fn path_edges(n: usize) -> Vec<(usize, usize)> {
        (1..n).map(|i| (i, i - 1)).collect()
    }

    fn binary_tree() -> LcaTree {
        // Complete binary tree with 7 nodes:
        //         0
        //        / \
        //       1   2
        //      / \ / \
        //     3  4 5  6
        LcaTree::new(7, vec![(1, 0), (2, 0), (3, 1), (4, 1), (5, 2), (6, 2)]).unwrap()
    }

    #[test]
    fn single_node() {
        let t = LcaTree::new(1, vec![]).unwrap();
        assert_eq!(t.lca(0, 0), 0);
        assert_eq!(t.parent(0), 0); // root is its own parent
    }

    #[test]
    fn path_graph_lca() {
        // Tree: 0 - 1 - 2 - 3 - 4 (path)
        let t = LcaTree::new(5, path_edges(5)).unwrap();

        // LCA on a path = the shallower node
        assert_eq!(t.lca(3, 4), 3);
        assert_eq!(t.lca(4, 3), 3);
        assert_eq!(t.lca(1, 4), 1);
        assert_eq!(t.lca(0, 4), 0);
        assert_eq!(t.lca(2, 2), 2);
        assert_eq!(t.lca(0, 0), 0);
    }

    #[test]
    fn star_graph_lca() {
        // Tree: 0 is root with children 1,2,3,4
        let edges = vec![(1, 0), (2, 0), (3, 0), (4, 0)];
        let t = LcaTree::new(5, edges).unwrap();

        assert_eq!(t.lca(1, 2), 0);
        assert_eq!(t.lca(3, 4), 0);
        assert_eq!(t.lca(0, 3), 0);
        assert_eq!(t.lca(2, 2), 2);
    }

    #[test]
    fn binary_tree_lca() {
        let t = binary_tree();

        assert_eq!(t.lca(3, 4), 1);
        assert_eq!(t.lca(5, 6), 2);
        assert_eq!(t.lca(3, 5), 0);
        assert_eq!(t.lca(4, 6), 0);
        assert_eq!(t.lca(3, 1), 1);
        assert_eq!(t.lca(1, 2), 0);
        assert_eq!(t.lca(0, 6), 0);
        assert_eq!(t.lca(3, 3), 3);

        assert_eq!(t.parent(3), 1);
        assert_eq!(t.parent(0), 0); // root is its own parent
    }

    #[test]
    fn lca_is_symmetric() {
        let t = binary_tree();
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t.lca(i, j), t.lca(j, i), "LCA({i},{j}) != LCA({j},{i})");
            }
        }
    }

    #[test]
    fn lca_with_self_is_self() {
        let t = binary_tree();
        for i in 0..7 {
            assert_eq!(t.lca(i, i), i);
        }
    }

    #[test]
    fn serialize_roundtrip() {
        let t = binary_tree();
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaTree::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n(), 7);
        for i in 0..7 {
            for j in 0..7 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
            }
            assert_eq!(t2.parent(i), t.parent(i));
        }
    }

    #[test]
    fn serialize_roundtrip_path() {
        let t = LcaTree::new(10, path_edges(10)).unwrap();
        let mut buf = Vec::new();
        t.serialize(&mut buf).unwrap();
        let t2 = LcaTree::load(&mut buf.as_slice()).unwrap();

        assert_eq!(t2.n(), 10);
        for i in 0..10 {
            for j in 0..10 {
                assert_eq!(t2.lca(i, j), t.lca(i, j));
            }
        }
    }

    // --- Validation error cases ---

    #[test]
    fn error_wrong_edge_count() {
        let err = LcaTree::new(4, vec![(1, 0), (2, 0)]).unwrap_err();
        assert!(err.contains("edges"), "{err}");
    }

    #[test]
    fn error_self_loop() {
        let err = LcaTree::new(3, vec![(1, 0), (1, 1)]).unwrap_err();
        assert!(err.contains("Self-loop") || err.contains("multiple parents"), "{err}");
    }

    #[test]
    fn error_multiple_parents() {
        // Node 2 gets two parents
        let err = LcaTree::new(4, vec![(1, 0), (2, 0), (2, 1)]).unwrap_err();
        assert!(err.contains("multiple parents"), "{err}");
    }

    #[test]
    fn error_disconnected() {
        // Nodes 1,2,3 form a cycle (each with one parent), node 0 has no children.
        // Edge count is correct (3 edges for 4 nodes), but 0 is unreachable as ancestor.
        let err = LcaTree::new(4, vec![(1, 2), (2, 3), (3, 1)]).unwrap_err();
        assert!(err.contains("disconnected") || err.contains("Root"), "{err}");
    }

    #[test]
    fn error_root_has_parent() {
        let err = LcaTree::new(3, vec![(1, 0), (0, 2)]).unwrap_err();
        assert!(err.contains("Root") || err.contains("parent"), "{err}");
    }

    #[test]
    fn error_out_of_range_node() {
        let err = LcaTree::new(3, vec![(1, 0), (5, 1)]).unwrap_err();
        assert!(err.contains("out-of-range"), "{err}");
    }

    #[test]
    fn error_zero_nodes() {
        let err = LcaTree::new(0, vec![]).unwrap_err();
        assert!(err.contains("at least one"), "{err}");
    }
}
