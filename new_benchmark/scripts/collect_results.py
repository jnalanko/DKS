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

    # In case of Bifrost I've added prints like this:
    # After simplify: 3.59465 seconds, 95416320 bytes RSS
    # After buildColors: 3.63468 seconds, 96509952 bytes RSS

    bifrost_time_after_simplify = None
    bifrost_time_after_write = None

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
                secs = float(s.split(":")[1])
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
        elif "After simplify:" in line:
            bifrost_time_after_simplify = float(line.split()[2])
        elif "After write:" in line:
            bifrost_time_after_write = float(line.split()[2])

    bifrost_coloring_time = None
    if bifrost_time_after_simplify != None and bifrost_time_after_write != None:
        bifrost_coloring_time = bifrost_time_after_write - bifrost_time_after_simplify
    return {"elapsed_seconds": elapsed_seconds, "max_rss_bytes": max_rss_bytes, "temp_disk": disk_peak, "bifrost_coloring_time": bifrost_coloring_time}

k_values = [15,23,31,39,47,55,63]

# Print tsv for plotting in R
print("\t".join(["k", "time_seconds", "mem_bytes"]))
for k in k_values:
    filename = f"logs/CHM13-k{k}-t1.log"
    try:
        res = parse_time_output(open(filename).readlines())
        time = res["elapsed_seconds"]
        print("{}\t{}\t{}".format(k, time, res["max_rss_bytes"]))
    except:
        pass



