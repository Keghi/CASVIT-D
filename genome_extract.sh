#!/bin/bash -x
#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com
while IFS= read -r line; do
 spec_nm=$(echo $line | cut -d ' ' -f1)
 genbank_acc=$(echo $line | cut -d ' ' -f2)
 unzip "../work/${genbank_acc}_crude.zip"
 mv "ncbi_dataset/data/${genbank_acc}" .
 rm -rf ncbi_dataset README.md
 python genomeheader_edit_2.py ${genbank_acc}
 #mv ${genbank_acc}/*.fna.modified ${genbank_acc}/*.fna
 mkdir -p ${spec_nm}
 /home/keshav/ucsc_tools/faToTwoBit ${genbank_acc}/*.fna.modified ${spec_nm}/${spec_nm}.2bit
 /home/keshav/ucsc_tools/twoBitInfo ${spec_nm}/${spec_nm}.2bit ${spec_nm}/chrom.sizes
 rm -fr ${genbank_acc} 
done < bla.txt
