use std::{cmp::{max, Reverse}, io::Write, ops::Range};
use jseqio::seq_db::SeqDB;
use crate::single_colored_kmers::{ColorVecValue, SingleColoredKmers};


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
    result: Vec<ColorVecValue>,

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
        let total_query_kmers = self.seqs.iter().fold(0_usize, |acc, rec| 
            acc + kmers_in_n(self.k, rec.seq.len()) 
        );
        let mut color_ids = Vec::<ColorVecValue>::with_capacity(total_query_kmers);

        for rec in self.seqs.iter() {
            for color in index.lookup_kmers(rec.seq, self.k) {
                color_ids.push(color);
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
    run_color: ColorVecValue, 
}

fn print_run<W: Write>(out: &mut W, seq_id: isize, run_color: ColorVecValue, range: Range<usize>) {
    // This code is almost duplicated in single_threaded_queries.rs
    if !range.is_empty() {
        match run_color {
            ColorVecValue::Single(c) => writeln!(out, "{seq_id}\t{}\t{}\t{}", range.start, range.end-1, c).unwrap(),
            ColorVecValue::Multiple => writeln!(out, "{seq_id}\t{}\t{}\t*", range.start, range.end-1).unwrap(),
            ColorVecValue::None => (), // Do not print runs of misses
        }
    }
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
                print_run(out, *cur_seq_id, *run_color, *p..(*p + *run_len));
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
                    print_run(out, *cur_seq_id, *run_color, *p .. (*p + *run_len));
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
        run_open: None, // Run of the same color (None counts as color). Value of None means that no run is active, not even a run of Nones.
        run_len: 0,
        run_color: ColorVecValue::None,
    };
    
    writeln!(out, "seq_rank\tfrom_kmer\tto_kmer\tcolor").unwrap();
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
        print_run(out, output_state.cur_seq_id, output_state.run_color, p..(p+output_state.run_len));
    }

    n_kmers_processed
    // Channel is dropped (= closed) here.
}

// Batch size is in nucleotides (= bytes)
pub fn lookup_parallel(n_threads: usize, mut queries: impl sbwt::SeqStream + Send, index: &SingleColoredKmers, batch_size: usize, k: usize, mut out: impl Write + Send) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(2); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<ProcessedQueryBatch>(2);

    let query_start = std::time::Instant::now();
    let mut n_bases_processed = 0_usize;

    let n_kmers_processed = std::thread::scope(|s| {

        let reader_handle = s.spawn(|| {

            // Reader thread that pushes batches for workers

            let mut batch = QueryBatch::new(0, k); // Initialize an empty batch
            while let Some(seq) = queries.stream_next() {
                let n = seq.len();
                let b = batch_size;

                crate::util::process_kmers_in_pieces(seq, k, batch_size, |piece_idx, piece: &[u8]|{
                    batch.push(piece, piece_idx > 0);

                    if batch.len() >= b {
                        // Swap the current batch with an empty batch, and send it to processing
                        let next_batch_id = batch.batch_id + 1;
                        let mut batch_to_send = QueryBatch::new(next_batch_id, k); // Empty batch
                        std::mem::swap(&mut batch, &mut batch_to_send);
                        batch_send.send(batch_to_send).unwrap();
                    }
                });

                n_bases_processed += n;
            }

            // Push the last remaining non-full batch. Can be empty but that's ok.
            batch_send.send(batch).unwrap();

            log::info!("All input read");
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
    
    log::info!("{} k-mers in {} base pairs queried in {} seconds (excluding index loading time)", n_kmers_processed, n_bases_processed, query_duration.as_secs());
    log::info!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / n_kmers_processed as f64);
    log::info!("Query time per base pair: {} nanoseconds", query_duration.as_nanos() as f64 / n_bases_processed as f64);

}

#[cfg(test)]
mod tests {
    use rand_chacha::rand_core::{RngCore, SeedableRng};
    use sbwt::{BitPackedKmerSortingMem, SeqStream};

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

    #[derive(Clone, Copy, Debug, Eq, PartialEq)]
    enum Color {
        Single(usize),
        Multiple,
    }

    fn random_test(batch_size: usize) {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(125);

        let mut sequences: Vec<Vec<u8>> = Vec::new();

        let k = 7;

        // Generate DNA sequences of random length. 
        // Each k-mer in the sequences must be unique! The following algorithm ensures that,
        // but it may loop forever if we get unlucky. The RNG seed is chosen so that this does not happen.
        let mut kmer_to_color = std::collections::HashMap::<Vec<u8>, Color>::new();
        for color in 0..2 {
            let len = (rng.next_u64() % 5000 + 1) as usize;
            let mut seq: Vec<u8> = Vec::new();
            while seq.len() < len {
                let base = rng.next_u64() % 4;
                let nucleotide = match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
                seq.push(nucleotide);
                if seq.len() >= k {
                    let kmer = &seq[seq.len()-k..];
                    if let Some(old) = kmer_to_color.get(kmer) {
                        if let Color::Single(old_color) = old {
                            if *old_color != color {
                                kmer_to_color.insert(kmer.to_vec(), Color::Multiple);
                            } // Else the color matches the old color -> keep as is
                        } // Else we have multiple colors -> keep as is
                    } else {
                        // Nothing stored yet -> store single color
                        kmer_to_color.insert(kmer.to_vec(), Color::Single(color));
                    }
                }
            }
            sequences.push(seq);
        }

        // Build SBWT
        // Use in-memory construction
        eprintln!("Building SBWT...");
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
        eprintln!("Building SingleColoredKmers...");
        let sck = SingleColoredKmers::new(sbwt, lcs, seqstreams, 3);
        eprintln!("SingleColoredKmers built");

        // Generate 1000 random queries of lengths between 1 and 100
        let mut queries: Vec<Vec<u8>> = Vec::new();
        let n_queries = 1000;
        for query_id in 0..n_queries {
            let len = rng.next_u64() % 100 + 1;
            let mut query: Vec<u8> = (0..len).map(|_| {
                let base = rng.next_u64() % 4;
                match base {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                }
            }).collect();

            if query_id == n_queries - 1 {
                // Make sure the last query contains a known k-mer so that the color last run remains open
                query.extend_from_slice(kmer_to_color.keys().next().unwrap());
            }
            queries.push(query);
        }

        let out_vec = Vec::<u8>::new();
        let mut out = std::io::Cursor::new(out_vec);
        lookup_parallel(2, MultiSeqStream::new(queries.clone()), &sck, batch_size, k, &mut out);

        // Parse output tsv line by line
        let output_str = String::from_utf8(out.into_inner()).unwrap();
        let output_lines = output_str.lines();
        // For each query, the starting positions and colors of found k-mers
        let mut found_kmers: Vec::<Vec::<(usize,Color)>> = vec![Vec::new(); queries.len()]; 
        for (line_idx, line) in output_lines.enumerate() {
            if line_idx == 0 { // tsv header
                assert_eq!(line, "seq_rank\tfrom_kmer\tto_kmer\tcolor");
            } else {
                let mut fields = line.split('\t');
                let seq_id: usize = fields.next().unwrap().parse().unwrap();
                let start: usize = fields.next().unwrap().parse().unwrap();
                let end: usize = fields.next().unwrap().parse().unwrap();
                let color_token = fields.next().unwrap();
                let color = if color_token == "*" {
                    Color::Multiple
                } else {
                    Color::Single(color_token.parse::<usize>().unwrap())
                };
                for i in start..=end {
                    found_kmers[seq_id].push((i, color));
                }
            }
        }

        // Verify against known answers
        for query_id in 0..queries.len() {
            let mut true_answer = Vec::<(usize, Color)>::new();
            for (i, kmer) in queries[query_id].windows(k).enumerate() {
                if let Some(color) = kmer_to_color.get(kmer) {
                    true_answer.push((i, *color));
                }
            }
            eprintln!("{:?}", true_answer);
            eprintln!("{:?}", found_kmers[query_id]);
            assert_eq!(true_answer, found_kmers[query_id]);
        }
    }

    #[test]
    fn run_random_tests() {
        random_test(7);
        random_test(8);
        random_test(1000);
    }
}