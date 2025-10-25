use std::{cmp::{max, Reverse}, io::{BufWriter, Write}, path::Path};

use jseqio::{reader::DynamicFastXReader, seq_db::SeqDB};

use crate::single_colored_kmers::SingleColoredKmers;


struct QueryBatch {
    seqs: SeqDB,

    batch_id: usize,
    sequence_starts: Vec<usize>, // Source sequence changes at these answer indices 
    chars_in_batch: usize,
    kmers_in_batch: usize,
    k: usize,
}

#[derive(Debug)]
struct ProcessedQueryBatch {
    result: Vec<Option<usize>>,

    batch_id: usize,
    sequence_starts: Vec<usize> 
}

impl QueryBatch {

    fn new(batch_id: usize, k: usize) -> Self {
        Self {
            seqs: SeqDB::new(),
            batch_id,
            sequence_starts: vec![],
            chars_in_batch: 0,
            kmers_in_batch: 0,
            k
        }
    }

    fn push(&mut self, seq: &[u8], extend_prev: bool) {
        if !extend_prev {
            self.sequence_starts.push(self.kmers_in_batch);
        }
        self.seqs.push_seq(seq);
        self.kmers_in_batch += kmers_in_n(self.k, seq.len());
        self.chars_in_batch += seq.len();

    }

    fn len(&self) -> usize {
        self.chars_in_batch
    }

    fn run(self, index: &SingleColoredKmers) -> ProcessedQueryBatch {
        let k = index.k();
        let total_query_kmers = self.seqs.iter().fold(0_usize, |acc, rec| 
            acc + kmers_in_n(k, rec.seq.len()) 
        );
        let mut color_ids = Vec::<Option::<usize>>::with_capacity(total_query_kmers);

        for rec in self.seqs.iter() {
            //eprintln!("{:?}", index.lookup_kmers(rec.seq).collect::<Vec::<Option::<usize>>>());
            for color in index.lookup_kmers(rec.seq) {
                if let Some(color) = color {
                    color_ids.push(Some(color));
                } else {
                    color_ids.push(None);
                }
            }
        }

        assert_eq!(color_ids.len(), total_query_kmers);
        ProcessedQueryBatch{
            result: color_ids,
            batch_id: self.batch_id,
            sequence_starts: self.sequence_starts,
        }
    }
}

impl PartialEq for ProcessedQueryBatch{
    fn eq(&self, other: &Self) -> bool {
        self.batch_id == other.batch_id
    }
}
impl Eq for ProcessedQueryBatch {}

impl PartialOrd for ProcessedQueryBatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other)) // Using the total order
    }
}

impl Ord for ProcessedQueryBatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.batch_id.cmp(&other.batch_id)
    }
}

fn kmers_in_n(k: usize, n: usize) -> usize {
    max(0, n as isize - k as isize + 1) as usize
}

struct OutputState {
    cur_seq_id: isize,
    run_open: Option<usize>,
    run_len: usize,
    run_color: Option<usize>, 
}

fn output_batch_result<W: Write>(batch: &ProcessedQueryBatch, state: &mut OutputState, out: &mut W) {

    let cur_seq_id = &mut state.cur_seq_id;
    let run_open = &mut state.run_open;
    let run_len = &mut state.run_len;
    let run_color = &mut state.run_color;

    let mut starts_ptr = batch.sequence_starts.iter().peekable();
    for (i, color) in batch.result.iter().enumerate() {
        while starts_ptr.peek().is_some_and(|&&s| s == i) {
            // New sequence starts. This closes the currently open run, if exists
            if let Some(p) = run_open {
                if let Some(c) = run_color {
                    // Write the run only if it's not None
                    writeln!(out, "{}\t{}\t{}\t{}", cur_seq_id, p, *p + *run_len - 1, c).unwrap();
                }
                *run_open = None;
                *run_len = 0;
            }
            starts_ptr.next();
            *cur_seq_id += 1;
        }

        match run_open {
            None => { 
                // We are at the start of a sequence -> open a new run
                *run_open = Some(0);
                *run_len = 1;
                *run_color = *color;
            },
            Some(p) => { 
                // See if we can extend the current run
                if *run_color == *color {
                    // Extend
                    *run_len += 1; 
                } else {
                    // Run ends
                    if let Some(c) = run_color {
                        // Write the run only if it's not None
                        writeln!(out, "{}\t{}\t{}\t{}", cur_seq_id, p, *p + *run_len-1, c).unwrap();
                    }
                    *run_open = Some(*p + *run_len);
                    *run_len = 1;
                    *run_color = *color;
                }
            }
        }
    }

}

// Returns the total number of k-mers in all the received batches
fn output_thread<W: Write>(query_results: crossbeam::channel::Receiver<ProcessedQueryBatch>, out: &mut W) -> usize {
    let mut batch_buffer = std::collections::BinaryHeap::<Reverse<ProcessedQueryBatch>>::new(); // Reverse makes it a min heap
    let mut n_kmers_processed = 0_usize;
    let mut next_batch_id = 0_usize;

    let mut output_state = OutputState {
        cur_seq_id: -1,
        run_open: None, // Run of the same color (None counts as color)
        run_len: 0,
        run_color: None,
    };
    
    while let Ok(batch) = query_results.recv() {
        batch_buffer.push(Reverse(batch)); // Reverse makes this a min heap

        loop { // Print all batches that can now be printed
            let min_batch = batch_buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.batch_id == next_batch_id {
                    output_batch_result(min_batch, &mut output_state, out);
                    n_kmers_processed += min_batch.result.len();
                    batch_buffer.pop();
                    next_batch_id += 1;
                } else {
                    break; // Not ready to print min_batch yet
                }
            } else {
                break; // Batch buffer is empty -> We are done
            }
        }

    }

    // All batches processed

    // The last run of the last batch remains open (unless it's closed by the start of a
    // sequence that has no k-mers). Let's write it.
    if let Some(p) = output_state.run_open {
        if let Some(c) = output_state.run_color {
            writeln!(out, "{}\t{}\t{}\t{}", output_state.cur_seq_id, p, p+output_state.run_len-1, c).unwrap();
        }
    }

    n_kmers_processed
    // Channel is dropped (= closed) here.
}

// Batch size is in nucleotides (= bytes)
pub fn lookup_parallel(n_threads: usize, mut queries: impl sbwt::SeqStream + Send, index: &SingleColoredKmers, batch_size: usize, mut out: impl Write + Send) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(2); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<ProcessedQueryBatch>(2);

    let query_start = std::time::Instant::now();

    let n_kmers_processed = std::thread::scope(|s| {

        let reader_handle = s.spawn(|| {

            // Reader thread that pushes batches for workers

            let mut batch = QueryBatch::new(0, index.k()); // Initialize an empty batch
            while let Some(seq) = queries.stream_next() {
                let n = seq.len();
                let b = batch_size;
                let k = index.k();

                if n < k {
                    // No full k-mers in this sequence -> push an empty sequence
                    batch.push(b"", false);
                    continue;
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

                    batch.push(piece, piece_idx > 0);

                    if batch.len() >= b {
                        let next_batch_id = batch.batch_id + 1;
                        batch_send.send(batch).unwrap();
                        batch = QueryBatch::new(next_batch_id, index.k());
                    }
                }
            }

            // Push the last remaining non-full batch. Can be empty but that's ok.
            batch_send.send(batch).unwrap();

            eprintln!("All input read");
            drop(batch_send); // Close the channel
        });

        let writer_handle = s.spawn(|| {
            let n_kmers = output_thread(output_recv, &mut out); // Returns number of k-mers processed
            out.flush().unwrap(); // They say this needs to be done because errors during drop are ignored
            n_kmers
        });

        let mut worker_handles = Vec::new();
        for _ in 0..n_threads {
            let output_send_clone = output_send.clone(); // Moved into worker
            let batch_recv_clone = batch_recv.clone(); // Moved into worker
            let index_ref = &index; // Moved into worker
            worker_handles.push(s.spawn(move || {
                while let Ok(batch) = batch_recv_clone.recv() {
                   let processed_batch = batch.run(index_ref);
                   output_send_clone.send(processed_batch).unwrap();
                }
                eprintln!("Worker done");
            }));
        }

        // Wait for threads to finish
        reader_handle.join().unwrap(); // All work batches pushed to workers
        worker_handles.into_iter().for_each(|w| w.join().unwrap()); // All batches processed
        drop(output_send); // All output written to channel -> can close the channel
        let n_kmers_processed = writer_handle.join().unwrap(); // All output written out

        #[allow(clippy::let_and_return)] // Is clearer to give it a name
        n_kmers_processed
    });

    let query_duration = query_start.elapsed();
    
    eprintln!("{} k-mers queried in {} seconds (excluding index loading time)", n_kmers_processed, query_duration.as_secs());
    eprintln!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / n_kmers_processed as f64);

}

#[cfg(test)]
mod tests {
    use sbwt::{BitPackedKmerSortingMem, SeqStream, StreamingIndex};

    use crate::{parallel_queries::lookup_parallel, single_colored_kmers::SingleColoredKmers};

    struct SingleSeqStream {
        seq: Vec<u8>,
        pos: usize,
    }

    impl SingleSeqStream {
        fn new(seq: Vec<u8>) -> Self {
            Self { seq, pos: 0 }
        }
    }

    impl SeqStream for SingleSeqStream {
        fn stream_next(&mut self) -> Option<&[u8]> {
            if self.pos == 0 {
                self.pos = 1;
                Some(&self.seq)
            } else {
                None
            }
        }
    }

    struct MultiSeqStream {
        seqs: Vec<Vec<u8>>,
        pos: usize,
    }

    impl MultiSeqStream {
        fn new(seqs: Vec<Vec<u8>>) -> Self {
            Self { seqs, pos: 0 }
        }
    }

    impl SeqStream for MultiSeqStream {
        fn stream_next(&mut self) -> Option<&[u8]> {
            if self.pos < self.seqs.len() {
                let ret = &self.seqs[self.pos];
                self.pos += 1;
                Some(ret)
            } else {
                None
            }
        }
    }

    #[test]
    fn build_testcase() {
        // Generate 1000 DNA sequences as byte slices, of random lengths between 1 and 100
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let mut sequences: Vec<Vec<u8>> = Vec::new();
        for _ in 0..1000 {
            let len = rng.gen_range(1..=100);
            let seq: Vec<u8> = (0..len).map(|_| {
                let base = rng.gen_range(0..4);
                match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            }).collect();
            sequences.push(seq);
        }

        // Build SBWT
        let k = 5;

        // Use in-memory construction
        let (sbwt, lcs) = sbwt::SbwtIndexBuilder::new()
            .add_rev_comp(false)
            .k(k)
            .build_lcs(true)
            .n_threads(3)
            .precalc_length(3)
            .algorithm(BitPackedKmerSortingMem::new().dedup_batches(false))
        .run_from_vecs(&sequences);
        let lcs = lcs.unwrap();

        let seqstreams: Vec<SingleSeqStream> = sequences.iter().map(|s| SingleSeqStream::new(s.clone())).collect();
        let sck = SingleColoredKmers::new(sbwt, lcs, seqstreams, 3);

        // Generate 1000 random queries of lengths between 1 and 100
        let mut queries: Vec<Vec<u8>> = Vec::new();
        for _ in 0..1000 {
            let len = rng.gen_range(1..=100);
            let query: Vec<u8> = (0..len).map(|_| {
                let base = rng.gen_range(0..4);
                match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            }).collect();
            queries.push(query);
        }

        let out_vec = Vec::<u8>::new();
        let mut out = std::io::Cursor::new(out_vec);
        lookup_parallel(2, MultiSeqStream::new(queries), &sck, 50, &mut out);

        eprintln!("Output:\n{}", String::from_utf8(out.into_inner()).unwrap());

    }
}