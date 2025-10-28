// Splits `seq` into pieces of length at most `max_piece_len`, such that the pieces overlap
// by k-1 characters. Calls f on pairs (piece_idx, piece). 
// If seq has length less than k, calls f with a single empty piece.
pub fn process_kmers_in_pieces<F: FnMut(usize, &[u8])>(seq: &[u8], k: usize, max_piece_len: usize, mut f: F) {
    let n = seq.len();
    let b = max_piece_len;

    if n < k {
        // No full k-mers in this sequence.
        f(0, b"");
        return;
    }

    // Let b be batch size and n be the length of the sequence.
    // Split the sequence into m pieces of length b except for the
    // last sequence that can have a shorter length, such that the
    // pieces overlap by k-1 characters and cover the whole sequence.
    // How many pieces will be there be? That is, what is the smallest
    // m such that b + (m-1)*(b-(k-1)) >= n? We must have b-k+1 > 0, or 
    // otherwise the inequality flips the wrong way around.
    // Assuming b-k+1 > 0, the solution is: m >= (n-k+1) / (b-k+1).
    // So we take the ceil of the righthand side.

    assert!(b as isize - k as isize + 1 > 0); // b-k+1 > 0
    let m = (n-k+1).div_ceil(b-k+1);
    assert!(m > 0); // This should be true since we checked for n < k earlier

    for piece_idx in 0..m {
        let pieces_before = piece_idx;
        let start = b*piece_idx - (k-1)*pieces_before;
        let piece = if piece_idx < m-1 {
            // Not the last piece: has full length b
            &seq[start..start+b]
        } else {
            &seq[start..] // Until the end (can have length shorter than b)
        };

        f(piece_idx, piece);
    }
}