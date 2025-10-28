#!/usr/bin/env python3
import sys
import re

def parse_time_output(lines):
    """
    Parse the output of /usr/bin/time -v.
    Returns a dict with 'elapsed_seconds' and 'max_rss_bytes'.
    """
    rss_re = re.compile(r"Maximum resident set size.*:\s*(\d+)")
    disk_re = re.compile(r"Temporary disk space peak: (\d+) bytes")

    elapsed_seconds = None
    max_rss_bytes = None
    disk_peak = 0

    for line in lines:
        # Match elapsed time — supports h:mm:ss, m:ss, or s.s formats
        if "Elapsed (wall clock) time" in line:
            # Two formats based on whether the time is more than a hour:
            # Elapsed (wall clock) time (h:mm:ss or m:ss): 1:13:01
            # Elapsed (wall clock) time (h:mm:ss or m:ss): 51:00.54
            s = line.split(" ")[-1].strip()
            if s.count(":") == 2: # h:mm:ss
                hours = int(s.split(":")[0])
                mins = int(s.split(":")[1])
                secs = int(s.split(":")[2])
                elapsed_seconds = hours * 60*60 + mins * 60 + secs
            else: # m:ss
                mins = int(s.split(":")[0])
                secs = float(s.split(":")[0])
                elapsed_seconds = mins * 60 + secs

        # Match max RSS (in kilobytes)
        elif "Maximum resident set size" in line:
            match = rss_re.search(line)
            if match:
                max_rss_bytes = int(match.group(1)) * 1024  # convert KB to bytes

        elif "Temporary disk space peak:" in line:
            match = disk_re.search(line)
            if match:
                disk_peak = int(match.group(1))

    return {"elapsed_seconds": elapsed_seconds, "max_rss_bytes": max_rss_bytes, "temp_disk": disk_peak}


def parse_lookup_stderr(lines):
    for line in lines:
        # Looking for a line like: 
        # Query time per base pair: 44.37868563833857 nanoseconds
        if "Query time per base pair" in line:
            return float(line.split(" ")[-2])
    return None

threads_list = [1,2,4,8,16,32]
k_list = [31,63]
mode_list = ["", "-mem"]
out = open("results.csv", 'w')
out.write("k,threads,time,mem,disk,mode,query_time_per_base_ns\n")

for t in threads_list:
    for k in k_list:
        for mode in mode_list:
            indexing_filename = f"logs/t{t}-k{k}{mode}.txt"
            query_filename = f"logs/query-t{t}-k{k}.txt"
            try:
                res = parse_time_output(open(indexing_filename).readlines())
                query_ns = parse_lookup_stderr(open(query_filename).readlines())
                res["query_time_per_base_ns"] = query_ns
                mode_string = "mem" if mode == "-mem" else "disk"
                out.write("{},{},{},{},{},{},{}\n".format(k,t,res["elapsed_seconds"],res["max_rss_bytes"],res["temp_disk"],mode_string, res["query_time_per_base_ns"]))
                print(f"Parsed {indexing_filename} {query_filename}")
            except:
                print(f"Warning: file {indexing_filename} or {query_filename} not found")
out.flush()
