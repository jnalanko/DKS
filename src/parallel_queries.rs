use std::{cmp::{max, min, Reverse}, io::Write, ops::Range, path::Path};

use jseqio::{reader::DynamicFastXReader, seq_db::SeqDB};

use crate::single_colored_kmers::SingleColoredKmers;


struct QueryBatch {
    seqs: SeqDB,

    batch_id: usize,
    sequence_starts: Vec<usize> // Source sequence changes at these answer indices 
}

#[derive(Debug)]
struct ProcessedQueryBatch {
    result: Vec<Option<usize>>,

    batch_id: usize,
    sequence_starts: Vec<usize> 
}

impl QueryBatch {
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

// Returns the total number of k-mers in all the received batches
fn output_thread<W: Write>(query_results: crossbeam::channel::Receiver<ProcessedQueryBatch>, out: &mut W) -> usize {
    let mut batch_buffer = std::collections::BinaryHeap::<Reverse<ProcessedQueryBatch>>::new(); // Reverse makes it a min heap
    let mut n_kmers_processed = 0_usize;
    let mut next_batch_id = 0_usize;
    let mut cur_seq_id = -1_isize;

    while let Ok(batch) = query_results.recv() {
        batch_buffer.push(Reverse(batch)); // Reverse makes this a min heap

        loop { // Print all batches that can now be printed
            let min_batch = batch_buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.batch_id == next_batch_id {
                    let mut starts_ptr = min_batch.sequence_starts.iter().peekable();
                    for (i, color) in min_batch.result.iter().enumerate() {
                        while starts_ptr.peek().is_some_and(|&&s| s == i) {
                            starts_ptr.next();
                            cur_seq_id += 1;
                            if cur_seq_id == 0 {
                                write!(out, "{cur_seq_id}").unwrap();
                            } else {
                                write!(out, "\n{cur_seq_id}").unwrap();
                            }
                        }

                        match color {
                            Some(color) => write!(out, " {color}").unwrap(),
                            None => write!(out, " X").unwrap(),
                        };
                    }

                    // In the last batch we can have sequence starts one past the end of the answers.
                    // In that case they are empty sequences. Print one line for each. 
                    for &s in starts_ptr {
                        assert!(s == min_batch.result.len());
                        cur_seq_id += 1;
                        if cur_seq_id == 0 {
                            write!(out, "{cur_seq_id}").unwrap();
                        } else {
                            write!(out, "\n{cur_seq_id}").unwrap();
                        }
                    }

                    n_kmers_processed += min_batch.result.len();
                    batch_buffer.pop();
                    next_batch_id += 1;
                } else {
                    break; // Not ready to print min_batch yet
                }
            } else {
                break; // Batch buffer is empty 
            }
        }
    }

    write!(out, "\n"); // One last newline to finish the last line

    //todo!(); // Need to add a newline after the very last seq
    n_kmers_processed
    // Channel is dropped (= closed) here.
}

// Batch size is in nucleotides (= bytes)
pub fn lookup_parallel(n_threads: usize, query_path: &Path, index: SingleColoredKmers, batch_size: usize) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(2); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<ProcessedQueryBatch>(2);

    let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
    let query_start = std::time::Instant::now();

    let n_kmers_processed = std::thread::scope(|s| {

        let reader_handle = s.spawn(|| {

            // Reader thread that pushes batches for workers

            let mut batch_seqs = SeqDB::new();
            let mut chars_in_batch = 0_usize;
            let mut kmers_in_batch = 0_usize;
            let mut seq_starts = Vec::<usize>::new();
            let mut batch_id = 0_usize;
            let mut n_seqs_read = 0_usize;

            while let Some(rec) = reader.read_next().unwrap() {
                n_seqs_read += 1;
                // Let b be batch size and n be the length of the sequence.
                // Split the sequence into m pieces of length b except for the
                // last sequence that can have a shorter length, such that the
                // pieces overlap by k-1 characters and cover the whole sequence.
                // Record a sequence break at the end of the last sequence.

                // How many pieces will be there be? That is, what is the smallest
                // m such that
                //
                // b + (m-1)*(b-(k-1)) >= n
                //
                // We must have b-k+1 > 0, or otherwise the inequality flips the wrong way around.
                // Assuming b-k+1 > 0, the solution is:
                //
                // m >= (n-k+1) / (b-k+1)
                //
                // So we take the ceil of the

                let seq = rec.seq; 
                let n = seq.len();
                let b = batch_size;
                let k = index.k();

                if n < k {
                    // No full k-mers in this sequence
                    if n_seqs_read > 1 {
                        seq_starts.push(kmers_in_batch);
                    }
                    continue;
                }

                assert!(b as isize - k as isize + 1 > 0); // b-k+1 > 0
                let m = (n-k+1).div_ceil(b-k+1);
                assert!(m > 0); // This should be true since we checked for n < k earlier

                // Sanity checks: 
                //     if n = k and b > k-1, then m = ceil(1 / positive), so m >= 1. Correct.
                //     if n = k+3 and b = k, then m = ceil(4) = 4. Correct.
                //     if n = 100, b = 30 and k = 10, then m = 5. Correct.
                // Seems to work.

                for piece_idx in 0..m {
                    let pieces_before = piece_idx;
                    let start = b*piece_idx - (k-1)*pieces_before;
                    let piece = if piece_idx < m-1 {
                        // Not the last piece: has full length b
                        &rec.seq[start..start+b]
                    } else {
                        &rec.seq[start..] // Until the end (can have length shorter than b)
                    };

                    if piece_idx == 0 {
                        seq_starts.push(kmers_in_batch);
                    }
                    batch_seqs.push_seq(piece);
                    kmers_in_batch += kmers_in_n(k, piece.len());
                    chars_in_batch += piece.len();

                    if chars_in_batch >= b {
                        let batch = QueryBatch {
                            seqs: batch_seqs,
                            batch_id,
                            sequence_starts: seq_starts,
                        };
                        batch_send.send(batch).unwrap();

                        // Reset the batch and tracking variables
                        batch_seqs = SeqDB::new();
                        seq_starts = Vec::new();
                        chars_in_batch = 0;
                        kmers_in_batch = 0;
                        batch_id += 1;
                    }
                }
            }

            eprintln!("Remaining seq starts: {:?}", seq_starts);

            // Push the last remaining non-full batch. Can be empty but that's ok.
            let batch = QueryBatch {
                seqs: batch_seqs,
                batch_id,
                sequence_starts: seq_starts,
            };
            batch_send.send(batch).unwrap();

            eprintln!("All input read");
            drop(batch_send); // Close the channel
        });

        let writer_handle = s.spawn(|| {
            let mut stdout = std::io::stdout();
            output_thread(output_recv, &mut stdout) // Returns number of k-mers processed
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
    
    /*for color in 0..index.n_colors() {
        let hits = color_hit_counts[color];
        println!("Color {}: {} hits ({:.2}%)", color, hits, hits as f64 / total_kmers_queried as f64 * 100.0);
    }
    */
    
    eprintln!("{} k-mers queried in {} seconds (excluding index loading time)", n_kmers_processed, query_duration.as_secs());
    //eprintln!("{:.2}% of query k-mers found", color_hit_counts.iter().sum::<usize>() as f64 / total_kmers_queried as f64 * 100.0);
    eprintln!("Query time per k-mer: {} nanoseconds", query_duration.as_nanos() as f64 / n_kmers_processed as f64);

}
