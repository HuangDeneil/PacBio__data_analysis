#!/bin/bash 

if [ ! -d "CARD_homolog_Data" ] ;then
    mkdir CARD_homolog_Data
fi
if [ ! -d "prevalence" ] ;then
    mkdir prevalence
fi

###  homolog fasta data
wget https://card.mcmaster.ca/latest/data
mv data CARD_homolog_Data/homolog_model.tar.bz2

###  all related data
wget https://card.mcmaster.ca/latest/variants
mv variants prevalence/prevalence.tar.bz2

## 解壓縮
cd CARD_homolog_Data
tar jxvf CARD_homolog_Data/homolog_model.tar.bz2
cd ..
cd prevalence
tar jxvf prevalence/prevalence.tar.bz2
gunzip *.gz
cd ..
