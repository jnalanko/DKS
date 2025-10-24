use std::path::Path;

use jseqio::reader::DynamicFastXReader;

use crate::single_colored_kmers::SingleColoredKmers;

pub fn lookup_single_threaded(query_path: &Path, index: &SingleColoredKmers){
    let mut reader = DynamicFastXReader::from_file(&query_path).unwrap();
    let mut seq_id = 0_usize;

    while let Some(rec) = reader.read_next().unwrap() {
        let mut run_start: Option<usize> = None;
        let mut run_color: Option<usize> = None;
        for (i, color) in index.lookup_kmers(rec.seq).enumerate() {
            if let Some(color) = color { // This k-mer has a color
                if let Some(s) = run_start {
                    // Currently on a run
                    if color == run_color.unwrap() {
                        // Ok. Extending the run
                    } else {
                        // Run ends -> print and open a new run
                        println!("{seq_id}\t{}\t{}\t{}", s, i-1, run_color.unwrap());
                        run_start = Some(i);
                        run_color = Some(color);
                    }
                } else {
                    // Currently not on a run -> open a run
                    run_start = Some(i);
                    run_color = Some(color);
                }
            } else { // This k-mer does not have a color
                // Terminate the current run, if exists
                if let Some(s) = run_start {
                    println!("{seq_id}\t{}\t{}\t{}", s, i-1, run_color.unwrap());
                    run_start = None;
                    run_color = None;
                }
            }
        }

        // Done processing the sequence. If there is still an open run, write it
        if let (Some(s), Some(c)) = (run_start, run_color) {
            println!("{seq_id}\t{}\t{}\t{}", s, rec.seq.len()-index.k(), c);
        } 

        seq_id += 1;
    }

}