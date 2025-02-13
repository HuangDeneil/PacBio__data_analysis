# [canu](https://github.com/marbl/canu)
## [document](https://canu.readthedocs.io/)
```bash
docker pull staphb/canu:2.2
docker run --rm -ti staphb/canu:2.2 canu -h

# Assembly
docker run --rm -ti \
-v $(pwd)/raw_data/:/raw_data/:ro \
-v $(pwd)/canu/:/data/ \
staphb/canu:2.2 \
canu genomeSize=3.2m \
 -p SAMN15455050 -d SAMN15455050 \
 -pacbio /raw_data/Pacbio/SRR12159828.fastq.gz

# Trim
docker run --rm -ti \
-v $(pwd)/raw_data/:/raw_data/:ro \
-v $(pwd)/canu/:/data/ \
staphb/canu:2.2 \
canu genomeSize=3.2m \
 -p SAMN15455050 -d SAMN15455050-erate-0.015 \
 correctedErrorRate=0.015 \
 -trimmed -corrected -pacbio SAMN15455050/SAMN15455050.trimmedReads.fasta.gz


```


# [quast](https://github.com/ablab/quast)
```bash
docker pull staphb/quast:5.3.0

docker run --rm -ti staphb/quast:5.3.0 quast.py


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

docker run --rm -ti staphb/bwa:0.7.18 bwa


## build index
time bwa index ${genome_fasta} -p $index_name

## Run alignment
time bwa mem -t 4 $index_name $fastq -o ${output}.sam


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

```bash
docker pull staphb/prokka:1.14.6
docker run --rm -ti staphb/prokka:1.14.6 prokka -h

docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/prokka:1.14.6 \
prokka \
 --outdir ./prokka/SAMN15455050 --force \
 --prefix SAMN15455050 --addgenes --locustag LAP \
 --increment 10 --gffver 2 --centre CDC --compliant \
 --genus Lactiplantibacillus --species plantarum \
 --kingdom Bacteria --gcode 11 --usegenus \
 --proteins /opt/prokka/db/trusted/Ecocyc-17.6 \
 --evalue 1e-9 --rfam \
 --cpus 4 \
 ./canu/canu.contigs.fasta

# Lactiplantibacillus plantarum


```


```bash
docker pull staphb/roary:3.13.0
docker run --rm -ti staphb/roary:3.13.0 roary -h



```

```bash
docker pull staphb/plasmidfinder:2.1.6_2024-03-07


```




```bash
docker pull staphb/circlator:1.5.5



```

## [igv-reports](https://github.com/igvteam/igv-reports)
```bash
docker pull staphb/igv-reports:1.12.0

```
