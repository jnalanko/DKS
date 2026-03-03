k_values = list(range(13,63+1,5))
n_threads = 6

for k in k_values:
    print(f"/usr/bin/time -v dks lookup -q data/CHM13.fasta.gz -i index/CHM13-k63.dks -t {n_threads} -k {k} 1> out/CHM13-k{k}.tsv 2> logs/CHM13-k{k}.log")
