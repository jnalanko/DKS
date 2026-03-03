# Data

The input is the T2T human genome assembly CHM13 version 2.0. Download from here: [https://github.com/marbl/CHM13](here). Currently, this [https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz](direct link) works. Put the file into `data/chm13v2.0.fa.gz`.

The file needs to be broken up into one file per sequence...

# Index construction

todo: external memory.

```
N_THREADS=4
dks build -k 63 -o index/CHM13-k63.dks -t $N_THREADS -i CHM13-fof.txt```
```
