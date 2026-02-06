use std::{ops::Range, path::Path};

use jseqio::reader::DynamicFastXReader;

use crate::single_colored_kmers::{ColorVecValue, SingleColoredKmers};

fn print_run(seq_id: usize, run_color: ColorVecValue, range: Range<usize>) {
    // This code is almost duplicated in parallel_queries.rs
    if !range.is_empty() {
        match run_color {
            ColorVecValue::Single(c) => println!("{seq_id}\t{}\t{}\t{}", range.start, range.end-1, c),
            ColorVecValue::Multiple => println!("{seq_id}\t{}\t{}\t*", range.start, range.end-1),
            ColorVecValue::None => (), // Do not print runs of misses
        }
    }
}

pub fn lookup_single_threaded(query_path: &Path, index: &SingleColoredKmers, k: usize){

    let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
    let mut seq_id = 0_usize;

    println!("seq_rank\tfrom_kmer\tto_kmer\tcolor");
    while let Some(rec) = reader.read_next().unwrap() {
        let mut run_start: usize = 0;
        let mut run_color: ColorVecValue = ColorVecValue::None;
        let mut n_kmers = 0;
        for (i, color) in index.lookup_kmers(rec.seq, k).enumerate() {
            if i == 0 { // Start a new run
                run_start = i;
                run_color = color;
            }
            else if color != run_color {
                // Run ends
                print_run(seq_id, run_color, run_start..i);
                run_start = i;
                run_color = color;
            } else {
                // Run continues
            }
            n_kmers += 1;
        }

        // Done processing the sequence. Close the last run 
        if n_kmers > 0 {
            print_run(seq_id, run_color, run_start..n_kmers);
        }

        seq_id += 1;
    }

}