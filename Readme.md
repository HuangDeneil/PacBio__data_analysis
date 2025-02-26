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

<!-- 
## [spades](https://github.com/ablab/spades)
```bash
docker pull staphb/spades:4.0.0


# A single paired-end library (separate files, gzipped):
bin/spades.py -1 left.fastq.gz -2 right.fastq.gz -o output_folder

``` 
-->

```bash
docker pull staphb/prokka:1.14.6
docker run --rm -ti staphb/prokka:1.14.6 prokka -h

docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/prokka:1.14.6 \
prokka \
 --outdir ./prokka/SAMN15455050 --force \
 --prefix SAMN15455050 --addgenes --locustag LAP \
 --increment 10 \
 --gffver 3 \
 --centre tig0000000 --compliant \
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



docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/roary:3.13.0 \
roary -e --mafft -p 4 roary/*.gff 

docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/roary:3.13.0 \
roary  roary/SAMN15455050.gff



docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/roary:3.12.0 \
roary  roary/SAMN15455050.gff

```

```bash
docker pull staphb/plasmidfinder:2.1.6_2024-03-07


```


<!-- 

```bash
docker pull staphb/circlator:1.5.5

docker run --rm -ti staphb/circlator:1.5.5 circlator plot -h 

circlator plot your_genome.fasta your_prokka_output.gff -o output_plot.pdf --colours CDS=blue,rRNA=green,tRNA=red --labels

docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/circlator:1.5.5 \
circlator plot \
./circlator/SAMN15455050.fasta \
./circlator/SAMN15455050.gff \
-o ./circlator/output_plot.pdf \
--colours CDS=blue,rRNA=green,tRNA=red --labels

``` -->


<!-- 
## [igv-reports](https://github.com/igvteam/igv-reports)
```bash
docker pull staphb/igv-reports:1.12.0
docker run --rm -ti staphb/igv-reports:1.12.0 create_report  -h 


docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/igv-reports:1.12.0 \
create_report 


docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/igv-reports:1.12.0 \
create_report \
--genome ./igv-report/SAMN15455050.fasta \
-a ./igv-report/SAMN15455050.gff \
-o ./igv-report/igv_report

docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/ \
staphb/igv-reports:1.12.0 \
igv-reports \
-g ./igv-report/SAMN15455050.fasta \
-a ./igv-report/SAMN15455050.gff \
-o ./igv-report/igv_report

igv-reports -g /data/your_genome.fasta -a /data/your_prokka_output.gff -o /data/igv_report

docker pull pegi3s/igv:2.13.0

docker run --rm -ti \
-e USERID=$UID \
-e USER=$USER \
-e DISPLAY=$DISPLAY \
-v ./igv/db:/var/db:Z \
-v ./tmp/.X11-unix:/tmp/.X11-unix \
--device /dev/dri/ \
-v "./igv/:/data" \
pegi3s/igv:2.13.0 sleep 1000000
``` 
-->


```bash
docker pull umigs/circleator:v1.0.2

docker run --rm -ti umigs/circleator:v1.0.2 circleator -h

docker run --rm -ti \
-v /home/share/SAMN15455050/:/home/circleator:rw \
umigs/circleator:v1.0.2 \
circleator \
 --data=prokka/SAMN15455050/SAMN15455050.gbk \
 --config=circleator/Circleator/conf/genes-percentGC-GCskew-1.cfg \
 > circleator/SAMN15455050.svg



docker run --rm -ti \
-v /home/share/SAMN15455050/:/home/circleator:rw \
umigs/circleator:v1.0.2 \
circleator \
 --config=circleator/Circleator/conf/genes-percentGC-GCskew-1.cfg \
 --data=prokka/SAMN15455050/SAMN15455050.gbk \

circleator
         --config=config-file.txt
         --data=annotation-and-or-sequence-data.gbk
        [--config_format=standard
         --data_dir=/some/place/to/look/for/files
         --sequence=eg-fasta-formatted-seq.fsa
         --seqlen=200000
         --contig_list=/path/to/tab-delim-contig-list.txt
         --contig_gap_size=20000
         --contig_min_size=5000
         --debug='all,coordinates,input,loops,misc,packing,tracks'
         --no_seq
         --rotate_degrees=0
         --scaled_segment_list='2000-3000:5,4000-5000:2,7000-9000:0.5'
         --scaled_segment_file=scaled-segments.txt
         --pad=400
         --log=/path/to/debug/logfile.txt
         --conf_dir=/etc/circleator/conf
         --help
         --man
         --version]

```

```bash
docker pull staphb/blast:2.16.0
docker run --rm -ti staphb/blast:2.16.0 makeblastdb -help
docker run --rm -ti staphb/blast:2.16.0 blastn -help


## build blastn DB
docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/:rw \
staphb/blast:2.16.0 \
makeblastdb \
 -dbtype nucl \
 -input_type fasta \
 -in CARD/canu.contigs.fasta \
 -title canu_assembly_genome \
 -parse_seqids \
 -out CARD/blastDB/contigs.blast.db



## blastn
docker run --rm -ti \
-v /home/share/SAMN15455050/:/data/:rw \
staphb/blast:2.16.0 \
blastn \
 -query ARGDB/2025Feb18/library/library.fna \
 -db CARD/blastDB/contigs.blast.db \
 -num_threads 4 \
 -out CARD/AMR.blastn.txt \
 -perc_identity 80 \
 -qcov_hsp_perc 80 \
 -evalue 1e-4 \
 -dust no \
 -soft_masking false  
 

# -outfmt 6
# -outfmt '6 bitscore nident pident qcovs evalue score sstart send' \
# -max_target_seqs 11
# -outfmt 17


```