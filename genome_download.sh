#!/bin/bash -P
#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com

while IFS= read -r line; do
 spec_nm=$(echo $line | cut -d ' ' -f1)
 genbank_acc=$(echo $line | cut -d ' ' -f2)
 wget -O "${genbank_acc}_crude.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${genbank_acc}/download?include_annotation_type=GENOME_FASTA,SEQUENCE_REPORT&filename=${genbank_acc}.zip" 
done < test1.txt



