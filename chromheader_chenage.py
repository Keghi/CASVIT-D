#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com

import os

with open('chromheader_change.txt', 'r') as f:
    lines = f.readlines()

with open('chromname_tbc.txt', 'w') as f:
    for line in lines:
        spec_name = line.split()[0]
        chrom = line.split()[1]
        with open(f'2bit_files/{spec_name}/chrom.sizes', 'r') as chrom_file:
            first_line = chrom_file.readline().strip()
        f.write(f'{spec_name} {chrom} {first_line}\n')
