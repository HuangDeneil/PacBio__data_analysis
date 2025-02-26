

nowPath=`pwd`
DBcodePath="/home/share/DB/05_ARGDB"
outputLibraryPath="/home/share/DB/ARGDB/2025Feb18/library"
outputBwaIndexPath="/home/share/DB/ARGDB/2025Feb18/bwa-db"
DBname="2025Feb18"


mkdir -p $outputLibraryPath
mkdir -p $outputBwaIndexPath

## download CARD database
bash $DBcodePath/download_CARD_data.sh

## make library.fna & ARO relation file
perl $DBcodePath/make_ARG_libraryFna.pl \
--inputFna $nowPath/CARD_homolog_Data/nucleotide_fasta_protein_homolog_model.fasta \
--ARO_index $nowPath/CARD_homolog_Data/aro_index.tsv \
--card_prevalence $nowPath/prevalence/card_prevalence.txt \
-o $outputLibraryPath




