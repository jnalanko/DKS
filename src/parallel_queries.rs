use std::{cmp::{max, Reverse}, io::Write, ops::Range, path::Path};

use jseqio::{reader::DynamicFastXReader, seq_db::SeqDB};

use crate::single_colored_kmers::SingleColoredKmers;


struct QueryBatch {
    seqs: SeqDB,

    batch_id: usize,
    sequence_breaks: Vec<usize> // List of nucleotide positions after which we should output a newline in the output
}

struct ProcessedQueryBatch {
    result: Vec<Option<usize>>,

    batch_id: usize,
    sequence_breaks: Vec<usize> // List of nucleotide positions after which we should output a newline in the output
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
            sequence_breaks: self.sequence_breaks,
        }
    }
}

impl ProcessedQueryBatch {
    fn write<W: Write>(&self, out: &mut W) {
        let mut seq_break_idx = 0;
        for (i, color) in self.result.iter().enumerate() {

            while seq_break_idx < self.sequence_breaks.len() && self.sequence_breaks[seq_break_idx] == i {
                out.write_all(b"\n").unwrap();
                seq_break_idx += 1;
            }

            match color {
                Some(color) => write!(out, "{color} ").unwrap(), // Todo: new space after last one
                None => write!(out, "X ").unwrap(), // Todo: new space after last one
            };
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

    while let Ok(batch) = query_results.recv() {
        batch_buffer.push(Reverse(batch)); // Reverse makes this a min heap

        loop { // Print all batches that can now be printed
            let min_batch = batch_buffer.peek();
            if let Some(min_batch) = min_batch {
                let min_batch = &min_batch.0; // Unwrap from Reverse
                if min_batch.batch_id == next_batch_id {
                    min_batch.write(out);
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

    n_kmers_processed
    // Channel is dropped (= closed) here.
}

// Batch size is in nucleotides (= bytes)
pub fn lookup_parallel(n_threads: usize, query_path: &Path, index: SingleColoredKmers, batch_size: usize) {
    let (batch_send, batch_recv) = crossbeam::channel::bounded::<QueryBatch>(2); // Read the next batch while the latest one is waiting to be processed
    let (output_send, output_recv) = crossbeam::channel::bounded::<ProcessedQueryBatch>(2);

    let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
    let query_start = std::time::Instant::now();
    let k = index.k();

    let n_kmers_processed = std::thread::scope(|s| {

        let reader_handle = s.spawn(|| {
            // Reader thread that pushes batches for workers
            let mut buf = SeqDB::new();
            let mut seq_breaks = Vec::<usize>::new();
            let mut n_chars_in_buf = 0_usize;
            let mut batch_id = 0_usize;
            while let Some(rec) = reader.read_next().unwrap() {
                if n_chars_in_buf + rec.seq.len() < batch_size {
                    buf.push_seq(rec.seq);
                    n_chars_in_buf += rec.seq.len();
                    assert!(n_chars_in_buf > 0);
                    seq_breaks.push(n_chars_in_buf-1);
                } else {
                    let mut head_len = batch_size - n_chars_in_buf;
                    let mut head = &rec.seq[0..head_len];
                    let mut tail = &rec.seq[head_len..];
                    loop { // Keep sending full batches for as long as we can

                        buf.push_seq(head);
                        n_chars_in_buf += head_len;

                        if tail.len() == 0 {
                            seq_breaks.push(n_chars_in_buf-1);
                        }

                        let batch = QueryBatch {
                            seqs: buf,
                            batch_id,
                            sequence_breaks: seq_breaks,
                        };

                        batch_send.send(batch);

                        batch_id += 1;
                        buf = SeqDB::new();
                        n_chars_in_buf = 0;
                        seq_breaks = Vec::new();

                        if tail.len() >= batch_size {
                            head = &tail[0..batch_size];
                            tail = &tail[batch_size..];
                            head_len = batch_size;
                        } else {
                            break;
                        }
                    }

                    if tail.len() > 0 {
                        buf.push_seq(tail);
                        n_chars_in_buf += tail.len();
                        seq_breaks.push(n_chars_in_buf-1);
                    }
                }
            }

            if n_chars_in_buf > 0 {
                // Last one
                seq_breaks.push(n_chars_in_buf - 1);
                let batch = QueryBatch {
                    seqs: buf,
                    batch_id,
                    sequence_breaks: seq_breaks,
                };
                batch_send.send(batch);
            }


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
