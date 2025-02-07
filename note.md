# [canu](https://github.com/marbl/canu)
## [document](https://canu.readthedocs.io/)
```bash
# Assemble
canu -correct \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -pacbio  pacbio.fastq

# trim
canu -trim \
  -p ecoli -d ecoli \
  genomeSize=4.8m \
  -corrected -pacbio ecoli/ecoli.correctedReads.fasta.gz

canu \
  -p ecoli -d ecoli-erate-0.039 \
  genomeSize=4.8m \
  correctedErrorRate=0.039 \
  -trimmed -corrected -pacbio ecoli/ecoli.trimmedReads.fasta.gz
```

```bash
docker pull staphb/canu:2.2
docker pull staphb/canu
docker run --rm -ti staphb/canu:2.2 canu -h


# Assembly
docker run --rm -ti \
-v $(pwd)/raw_data/:/raw_data/:ro \
-v $(pwd)/canu/:/data/ \
staphb/canu:2.2 \
canu \
  -p SAMN15455050 -d SAMN15455050 \
  genomeSize=3.2m \
  -pacbio /raw_data/Pacbio/SRR12159828.fastq.gz

# Trim
docker run --rm -ti \
-v $(pwd)/raw_data/:/raw_data/:ro \
-v $(pwd)/canu/:/data/ \
staphb/canu:2.2 \
canu \
 -p SAMN15455050 -d SAMN15455050-erate-0.0015 \
 genomeSize=3.2m \
 correctedErrorRate=0.039 \
 -trimmed -corrected -pacbio SAMN15455050/SAMN15455050.trimmedReads.fasta.gz


```


# [quast](https://github.com/ablab/quast)
```bash
docker pull staphb/quast:5.3.0

./quast.py test_data/contigs_1.fasta \
           test_data/contigs_2.fasta \
        -r test_data/reference.fasta.gz \
        -g test_data/genes.txt \
        -1 test_data/reads1.fastq.gz -2 test_data/reads2.fastq.gz \
        -o quast_test_output

```


## [busco](https://gitlab.com/ezlab/busco)
```bash
docker pull staphb/busco:5.8.2

```

## [bwa](https://github.com/lh3/bwa)
```bash
docker pull staphb/bwa:0.7.18

time bwa index i -p argv$output_path/sample_id


time bwa mem -t 4 $index $fastq -o ${output}.sam


```

## [samtools](https://github.com/samtools/samtools)
```bash
docker pull staphb/samtools:1.21

# extracting umap reads
time samtools fastq -@ 4 -f 4 ${i} > ${output_path}/sample_id.unmap.fq 

```


## [spades](https://github.com/ablab/spades)
```bash
docker pull staphb/spades:4.0.0


# A single paired-end library (separate files, gzipped):
bin/spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder

```



