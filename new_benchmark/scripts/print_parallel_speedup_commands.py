k = 63

for n_threads in [1,2,3,4,5,6,7,8]:
    print(f"/usr/bin/time -v dks lookup -q data/CHM13.fasta.gz -i index/CHM13-k63.dks -t {n_threads} -k {k} 1> out/CHM13-k{k}-t{n_threads}.tsv 2> logs/CHM13-k{k}-t{n_threads}.log")
