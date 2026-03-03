#!/usr/bin/env python3
import sys
import re

def parse_time_output(lines):

    elapsed_seconds = None
    query_time_ns_per_bp = None
    max_rss_bytes = None

    rss_re = re.compile(r"Maximum resident set size.*:\s*(\d+)")

    for line in lines:
        # Match elapsed time — supports h:mm:ss, m:ss, or s.s formats
        if "Elapsed (wall clock) time" in line:
            s = line.split(" ")[-1].strip()
            if s.count(":") == 2: # h:mm:ss
                hours = int(s.split(":")[0])
                mins = int(s.split(":")[1])
                secs = int(s.split(":")[2])
                elapsed_seconds = hours * 60*60 + mins * 60 + secs
            else: # m:ss
                mins = int(s.split(":")[0])
                secs = float(s.split(":")[1])
                elapsed_seconds = mins * 60 + secs

        # Match max RSS (in kilobytes)
        elif "Maximum resident set size" in line:
            match = rss_re.search(line)
            if match:
                max_rss_bytes = int(match.group(1)) * 1024  # convert KB to bytes
        elif "Query time per base pair" in line:
            query_time_ns_per_bp = float(line.split()[-2])

    return {"elapsed_seconds": elapsed_seconds, "max_rss_bytes": max_rss_bytes, "query_time_ns_per_bp": query_time_ns_per_bp}

k = 63
n_threads_values = [1,2,3,4,5,6,7,8]

# Print tsv for plotting in R
print("\t".join(["threads", "elapsed_seconds", "query_time_ns_per_bp", "mem_bytes"]))
for n_threads in n_threads_values:
    filename = f"logs/CHM13-k{k}-t{n_threads}.log"
    try:
        res = parse_time_output(open(filename).readlines())
        mem = res["max_rss_bytes"]
        time = res["elapsed_seconds"]
        ns_per_bp = res["query_time_ns_per_bp"]
        print("{}\t{}\t{}\t{}".format(n_threads, time, ns_per_bp, mem))
    except:
        sys.stderr.write("Error parsing " + filename + "\n")
