#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com

import os
import collections
import logging
import subprocess
import argparse
import logging

# set up logging

logger = logging.getLogger()
logger.setLevel(logging.INFO)

file_handler = logging.FileHandler('gene_checker_logfile.txt')
file_handler.setLevel(logging.INFO)

logger.addHandler(file_handler)

#funtion for input  
def parse_args():
    parser = argparse.ArgumentParser(description='Extract fasta files for the gene present in the file')
    parser.add_argument('input_file', help='Input file containing gene information')
    
    return parser.parse_args()

#function to Extract maf of desired locations from the bigger .maf files
def run_maf_extract(gene_name, chr_number, exon_coords):
    for exon_number, exon_coord in enumerate(exon_coords, start=1):
        maf_file = f"{gene_name}_exon{exon_number}.maf"
        if os.path.exists(maf_file):
            print(f"{maf_file} already exists, skipping mafExtract command")
        else:
            mafExtract_cmd = f"mafExtract multiz470way_maf/mutiz470.bb -region={chr_number}:{exon_coord} {maf_file}"
            process = subprocess.run(mafExtract_cmd, shell=True, check=True)  # sanity check
            if process.returncode == 0:
                print(f"mafExtract command for {gene_name} exon {exon_number} completed successfully and {maf_file} is created")
            else:
                print(f"mafExtract command for {gene_name} exon {exon_number} failed with return code {process.returncode}")
    
    return maf_file

#function to generate the cordinates of the exons in .maf files
def generate_sequence_dict(gene_name):
   
        sequence_data = collections.defaultdict(dict)
        maf_file = f'{gene_name}_exon{exon_number}.maf'

        #calculates the start and end of the block(MSAB)
        with open(maf_file, 'r') as f:
            lines = (line.strip() for line in f.readlines() if line.strip())
        for line in lines:
            if line.startswith('s'):
                fields = line.split()
                seq_name, start, size, strand, src_size, maf_seq = fields[1:]
                species_name, chr_number = seq_name.split('.', 1)
                if strand == '+':
                    start = int(start) 
                    end = int(start) + int(size) 
                elif strand == '-':
                    end = int(src_size) - int(start) 
                    start = end - int(size) 
                    
                if 'seq_start' not in sequence_data[species_name]:
                    sequence_data[species_name].update({
                        'chr_number': chr_number,
                        'strand': strand,
                        'seq_start': int(start),
                        'seq_end': int(end)
                    })
                elif int(start) < sequence_data[species_name]['seq_start']:
                    sequence_data[species_name]['seq_start'] = int(start)
                elif int(end) > sequence_data[species_name]['seq_end']:
                    sequence_data[species_name]['seq_end'] = int(end)  

        #creates the cordinate table of the gene for each speices 
        coordinate_table = f'{gene_name}_exon{exon_number}_coordinate_table.txt'
        with open(coordinate_table, 'w') as outfile:
            outfile.write('{:<15} {:<15} {:<15} {:<15} {:<15}\n'.format('species Name', 'Chromosome', 'Strand', 'Sequence Start', 'Sequence End'))
            for species_name, data in sequence_data.items():
                outfile.write('{:<15} {:<15} {:<15} {:<15} {:<15}\n'.format(species_name, data['chr_number'], data['strand'], data['seq_start'], data['seq_end']))
        print(f'{coordinate_table} created')
        
        return sequence_data

def create_species_cord_table(gene_name, exon_number_list):
    # create an empty dictionary to store the coordinates of each species for each exon
    species_cord = {}
    # loop through each exon and generate the coordinate table for each species
    for exon_number in exon_number_list:
        coordinate_table = f'{gene_name}_exon{exon_number}_coordinate_table.txt'
        # read the coordinate table file and store the data in the dictionary
        with open(coordinate_table, 'r') as f:
            next(f) # skip header line
            for line in f:
                species_name, chr_number, strand, seq_start, seq_end = line.strip().split()
                # if this is the first time seeing this species, create a new dictionary entry for it
                if species_name not in species_cord:
                    species_cord[species_name] = {'Chromosome': [], 'Strand': [], 'Sequence_Start': [], 'Sequence_End': []}
                # add the coordinates for this exon to the species's dictionary entry
                species_cord[species_name]['Chromosome'].append(chr_number)
                species_cord[species_name]['Strand'].append(strand)
                species_cord[species_name]['Sequence_Start'].append(seq_start)
                species_cord[species_name]['Sequence_End'].append(seq_end)

        # write the species coordinate table to a file
    species_cord_file = f'species_cord_{gene_name}.txt'
    with open(species_cord_file, 'w') as f:
        # write the header line
        #f.write('##species Name##\n')
        # write the exon names in the order given in exon_number_list
        #f.write('Exons\t')
        #f.write('\n')
        # write the coordinates for each species
        for species_name, data in species_cord.items():
            f.write(f'###{species_name}###\n')
            for key, values in data.items():
                #f.write(f'{key}\t')
                # write the coordinates in the order given in values
                for value in values:
                    f.write(f'{value}\t')
                f.write('\n')

    print(f'{species_cord_file} created')



def remove_files(gene_name, exon_numbers):
    for exon_number in exon_numbers:
        maf_file = f"{gene_name}_exon{exon_number}.maf"
        coordinate_table = f'{gene_name}_exon{exon_number}_coordinate_table.txt'
        os.remove(maf_file)
        os.remove(coordinate_table)

#The main script     
if __name__ == '__main__':
    args = parse_args()

    with open(args.input_file, 'r') as f:  
        for line in f:
            split_line = line.strip().split()
            gene_name, chr_number, gene_strand = split_line[:3]
            exon_coord_str = split_line[3:]
            exon_coords = exon_coord_str[0].split(',')
            num_exons = len(exon_coords)
            if gene_strand == '-':
                exon_coords.reverse()
                exon_numbers = list(range(num_exons, 0, -1))
            else:
                exon_numbers = list(range(1, num_exons+1))
            print(f"Processing gene {gene_name} on chromosome {chr_number} strand {gene_strand} with {num_exons} exons and coordinates {exon_coords}")

            run_maf_extract(gene_name, chr_number, exon_coords)

            sequence_data = {}
            for exon_number in exon_numbers:
                generate_sequence_dict(gene_name)
            
            create_species_cord_table(gene_name, exon_numbers)
            remove_files(gene_name, exon_numbers)
