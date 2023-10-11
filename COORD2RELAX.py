import os
import collections
import logging
import subprocess
import argparse
from Bio.Seq import Seq #for reverse compilmenting sequences
from Levenshtein import distance #match the sequence in the maf_gene_seqs to gene_seq
logging.basicConfig(filename='log_file.txt', level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')


#funtion for input  
def parse_args():
    parser = argparse.ArgumentParser(description='Extract fasta files for the gene present in the file')
    parser.add_argument('input_file', help='Input file containing gene information')
    
    return parser.parse_args()

'''
#funtion to create the index file
def run_maf_index(chr_number):
    maf_bb_file = f"chr{chr_number}_maf.bb"
    if not os.path.exists(maf_bb_file):
        mafIndex_cmd = f"mafIndex -chromSizes=chrom.sizes chr{chr_number}.maf {maf_bb_file}"
        subprocess.run(mafIndex_cmd, shell=True, check=True)
    else:
        print(f"{maf_bb_file} already exists, skipping mafIndex command")
        
    return maf_bb_file
'''

#function to Extract maf of desired locations from the bigger .maf files
def run_maf_extract(gene_name, chr_number, exon_coords):
    for exon_number, exon_coord in enumerate(exon_coords, start=1):
        maf_file = f"{gene_name}_exon{exon_number}.maf"
        if os.path.exists(maf_file):
            print(f"{maf_file} already exists, skipping mafExtract command")
        else:
            mafExtract_cmd = f"mafExtract  multiz470way_maf/mutiz470.bb -region=chr{chr_number}:{exon_coord} {maf_file}"
            print(mafExtract_cmd)
            process = subprocess.run(mafExtract_cmd, shell=True, check=True)  # sanity check
            if process.returncode == 0:
                print(f"mafExtract command for {gene_name} exon {exon_number} completed successfully and {maf_file} is created")
            else:
                print(f"mafExtract command for {gene_name} exon {exon_number} failed with return code {process.returncode}")
    
    return maf_file


#function to generate the cordinates of the exons in .maf files
def generate_coordinate_tables(gene_name):
        sequence_data = {}
        prev_rStatus = prev_rCount = None
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
                    
                #store sequence for each species 
                if species_name not in maf_gene_seqs:
                    maf_gene_seqs[species_name] = ""
                maf_gene_seqs[species_name] += maf_seq.replace('\n', '').replace('-', '').upper()

                #store strand details of each species
                if species_name not in exon_strands:
                    exon_strands[species_name] = {}
                exon_strands[species_name][f"exon{exon_number}"] =  strand
                    
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
            
            #count the missing bases in the maf file
            elif line.startswith('i'):
                fields = line.split()
                seq_name, lStatus, lCount, rStatus, rCount = fields[1:]
                species_name, chr_number = seq_name.split('.', 1)
                if species_name not in missing_bases_maf:
                    missing_bases_maf[species_name] = 0
                if lStatus == "I" and prev_rStatus == lStatus and prev_rCount == lCount:
                    missing_bases_maf[species_name] += int(lCount)
                prev_rStatus = rStatus
                prev_rCount = rCount
                

        #creates the cordinate table of the gene for each speices 
        coordinate_table = f'{gene_name}_exon{exon_number}_coordinate_table.txt'
        with open(coordinate_table, 'w') as outfile:
            outfile.write('{:<15} {:<15} {:<15} {:<15} {:<15}\n'.format('species Name', 'Chromosome', 'Strand', 'Sequence Start', 'Sequence End'))
            for species_name, data in sequence_data.items():
                outfile.write('{:<15} {:<15} {:<15} {:<15} {:<15}\n'.format(species_name, data['chr_number'], data['strand'], data['seq_start'], data['seq_end']))
        print(f'{coordinate_table} created')
        
        return sequence_data

'''
# fucntion to create fasta sequnuces of exons       
def generate_sequence_dict(sequence_data, gene_name, exon_number):
    for species_name, species_data in sequence_data.items():
        spec_exon_fa = f'{gene_name}_{species_name}_exon{exon_number}.fa'
        command = f'twoBitToFa 2bit_files/{species_name}/{species_name}.2bit:{species_data["chr_number"]}:{species_data["seq_start"]}-{species_data["seq_end"]} {spec_exon_fa}'
        process = subprocess.run(command, shell=True, check=True)
        if process.returncode == 0:
            print(f"twoBitToFa command for {gene_name} {species_name} exon {exon_number} completed successfully and {spec_exon_fa} is created")
        else:
            print(f"twoBitToFa command for {gene_name} {species_name} exon {exon_number} failed with return code {process.returncode}")
                
        with open(spec_exon_fa, 'r') as f:
            lines = f.readlines()
            seq = Seq(''.join(lines[1:]).strip())
            if species_data["strand"] == '-':
                rc_seq = seq.reverse_complement()
                gene_seq = str(rc_seq).replace('\n', '').upper()
            else:
                gene_seq = str(seq).replace('\n', '').upper()
        if species_name not in gene_seqs:
            gene_seqs[species_name] = gene_seq
        else:
            gene_seqs[species_name] += gene_seq 
            
        if species_name == list(sequence_data.keys())[0]:
            if exon_number not in query_seq:
                query_seq[exon_number] = {}
            query_seq[exon_number][f'>{gene_name}_{species_name}_exon{exon_number}'] = gene_seq

'''

# fucntion to create fasta sequnuces of exons       
def generate_sequence_dict(sequence_data, gene_name, exon_number):
    for species_name, species_data in sequence_data.items():
        spec_exon_fa = f'{gene_name}_{species_name}_exon{exon_number}.fa'
        command = f'twoBitToFa 2bit_files/{species_name}/{species_name}.2bit:{species_data["chr_number"]}:{species_data["seq_start"]}-{species_data["seq_end"]} {spec_exon_fa}'
        try:
            process = subprocess.run(command, shell=True, check=True)
            if process.returncode == 0:
                print(f"twoBitToFa command for {gene_name} {species_name} exon {exon_number} completed successfully and {spec_exon_fa} is created")
            else:
                print(f"twoBitToFa command for {gene_name} {species_name} exon {exon_number} failed with return code {process.returncode}")
        except subprocess.CalledProcessError as e:
            print(f"twoBitToFa command for {gene_name} {species_name} exon {exon_number} failed: {str(e)}")
                
        with open(spec_exon_fa, 'r') as f:
            lines = f.readlines()
            seq = Seq(''.join(lines[1:]).strip())
            if species_data["strand"] == '-':
                rc_seq = seq.reverse_complement()
                gene_seq = str(rc_seq).replace('\n', '').upper()
            else:
                gene_seq = str(seq).replace('\n', '').upper()
        if species_name not in gene_seqs:
            gene_seqs[species_name] = gene_seq
        else:
            gene_seqs[species_name] += gene_seq 
            
        if species_name == list(sequence_data.keys())[0]:
            if exon_number not in query_seq:
                query_seq[exon_number] = {}
            query_seq[exon_number][f'>{gene_name}_{species_name}_exon{exon_number}'] = gene_seq

    return species_name

#function to check wheather the sequences obtained from 2bit files are identical to sequnces in maf            
def check_sequence_similarity(gene_seqs, maf_gene_seqs, threshold=0.95):
    all_identical = True
    for species, seq in gene_seqs.items():
        if species in maf_gene_seqs:
            maf_seq = maf_gene_seqs[species]
            dist = distance(seq, maf_seq)
            similarity_ratio = 1 - dist / max(len(seq), len(maf_seq))
            if similarity_ratio >= threshold:
                pass
            else:
                all_identical = False
                print(f"{species} sequence is not identical to sequence in maf file")
        else:
            print(f"No maf_gene_seq found for {species}")    
    if all_identical:
        print("All sequences are identical to sequence in maf file")
        

#function to create the mutiple sequence file and add the query sequnce as required by cesar software
def write_ms_file(cesar_input_file, query_seq):
    with open(cesar_input_file, 'w') as ms_file:
        for exon_number in sorted(query_seq.keys()):
            for header, seq in query_seq[exon_number].items():
                if gene_strand == '-':
                    seq = Seq(seq).reverse_complement()
                if exon_number == 1:
                    ms_file.write(f'{header}_exon{exon_number}\textra/tables/human/firstCodon_profile.txt\textra/tables/human/do_profile.txt\n')
                elif exon_number == len(query_seq):
                    ms_file.write(f'{header}_exon{exon_number}\textra/tables/human/acc_profile.txt\textra/tables/human/lastCodon_profile.txt\n')
                else:
                    ms_file.write(f'{header}_exon{exon_number}\textra/tables/human/acc_profile.txt\textra/tables/human/do_profile.txt\n')
                
                ms_file.write(f'{seq}\n')
    print(f'{cesar_input_file} is created and {gene_name}_hg38 sequence is added')

#comapre gene lenghts from 2bit files and maf files to know weather they are accurate or not
def compare_gene_lengths(gene_seqs, maf_gene_seqs, missing_bases_maf):
    for species_name in gene_seqs:
        if species_name not in maf_gene_seqs or species_name not in missing_bases_maf:
            continue
        if species_name in maf_gene_seqs:
            gene_len = len(gene_seqs[species_name])
            maf_gene_len = len(maf_gene_seqs[species_name].replace('\n','')) + int(missing_bases_maf[species_name])
            if gene_len == maf_gene_len:
                print(f"{species_name}: Gene length ({gene_len}) matches with MAF gene length ({maf_gene_len})")
            else:
                print(f"{species_name}: Gene length ({gene_len}) does not match with MAF gene length ({maf_gene_len})")
        

#function to write fasta files from the gene_seqs dictionary
def write_fasta_files(gene_name, gene_strand, gene_seqs, cesar_input_file, exon_strands):
    with open(cesar_input_file, "r") as f:
        result = []
        start_letter = 0
        for line in f:
            if line.startswith('>'):
                result.append(line)
            else:
                seq = line.rstrip()
                seq_len = len(seq)-start_letter
                rem = seq_len % 3
                #this to make the first letters small
                if start_letter == 0 or start_letter == 3:
                    pass
                elif start_letter == 1:
                    seq = seq[:1].lower() + seq[1:]
                elif start_letter == 2:
                    seq = seq[:2].lower() + seq[2:]
                #this to make the last letters small
                if rem == 1:
                    seq = seq[:-1] + seq[-1].lower()
                elif rem == 2:
                    seq = seq[:-2] + seq[-2:].lower()
                result.append(seq + '\n')
                start_letter = 3 - rem  
            
        result.append("####" + '\n' )

    #appends the remaining sequnces to result    
    for species_name, gene_seq in gene_seqs.items():
        if species_name != list(gene_seqs.keys())[0]:
            if species_name != list(gene_seqs.keys())[0]:
                strand_values = exon_strands.get(species_name, {}).values()
                if '-' in strand_values and '+' in strand_values:
                    print(f"Skipping {species_name} because not all exon strands have the same strand.")
                    logging.info(f"Skipping {gene_name}_{species_name} because not all exon strands have the same strand.")
                    continue
                else:
                    header = f'>{gene_name}_{species_name}\n'
                    if gene_strand == '-':
                        seq = Seq(gene_seq)
                        rc_seq = seq.reverse_complement()
                        gene_seq = str(rc_seq)

                    if gene_seq:  # Check if gene_seq is not empty
                        result.append(header + gene_seq + '\n')
                    else:
                        logging.info(f"Skipping {gene_name}_{species_name} because of an empty sequence.")
                        
    #write the updated sequences to the file
    with open(cesar_input_file, "w") as f:
        f.writelines(result)
        print(f'Other species fasta sequences are updated in the {cesar_input_file}')   
'''  
    #appends the remaining sequnces to result 
    for species_name, gene_seq in gene_seqs.items():
        if species_name != list(gene_seqs.keys())[0]:
            if species_name != list(gene_seqs.keys())[0]:
                strand_values = exon_strands.get(species_name, {}).values()
                if '-' in strand_values and '+' in strand_values:
                    print(f"Skipping {species_name} because not all exon strands have the same strand.")
                    logging.info(f"Skipping {gene_name}_{species_name} because not all exon strands have the same strand.")
                    continue
                else:
                    header = f'>{gene_name}_{species_name}\n'
                    if gene_strand == '-':
                        seq = Seq(gene_seq)
                        rc_seq = seq.reverse_complement()
                        gene_seq = str(rc_seq)
                        result.append(header + gene_seq + '\n')
                    else:
                        result.append(header + gene_seq + '\n')

    #write the updated sequences to the file
    with open(cesar_input_file, "w") as f:
        f.writelines(result)
        print(f'Other species fasta sequences are updated in the {cesar_input_file}')
'''     

#function to run ceasar and create MSA file       
def create_MSA_file(cesar_input_file, gene_name):
    command = f'/home/SharmaLab/CESAR2.0/cesar {cesar_input_file} > {cesar_out_file}'
    process = subprocess.run(command, shell=True, check=True)
    if process.returncode == 0:
        print(f"./cesar command for {gene_name} completed successfully and {cesar_out_file} is created")
    else:
        print(f"./cesar command for {gene_name} failed with return code {process.returncode}")



#function to covert cesar output to multiple sequece alignment
def process_cesar_out_file(cesar_out_file, gene_name):
    with open(cesar_out_file, 'r') as f:
        lines = f.readlines()
        result = []
        for i in range(0, len(lines), 4):
            ref_header = lines[i].strip()
            ref_seq = lines[i+1].replace(' ', '-').strip()
            query_header = lines[i+2].strip()
            query_seq = lines[i+3].replace(' ', '-').strip()

            new_ref_seq = ''
            new_query_seq = ''
            for ref_char, query_char in zip(ref_seq, query_seq):
                if ref_char == '>':
                    continue
                elif ref_char == ' ':
                    new_ref_seq += '-'
                    new_query_seq += query_char
                else:
                    new_ref_seq += ref_char if ref_char != ' ' else '-'
                    new_query_seq += query_char

            if len(new_ref_seq) % 3 != 0:
                continue  # Skip the block if the length is not a multiple of 3

            species_name = query_header.split('_')[1]
            padding = len(species_name) - len('hg38')
            padding_str = ' ' * padding

            maf_header = f'##maf version=1 scoring=tba.v8'
            maf_ref_seq = f's hg38.chrX{padding_str} 1 1 + 1 {new_ref_seq}'
            maf_query_seq = f's {species_name}.chrX 1 1 + 1 {new_query_seq}'

            result.append(maf_header)
            result.append(maf_ref_seq)
            result.append(maf_query_seq)

        result = [line for line in result if line.count('N') <= 30 and line.count('n') <= 30]

        # Write the output to a file
        maf_join_input = f'{gene_name}_maf_join_input.maf'

        with open(maf_join_input, 'w') as output_file:
            output_file.write('\n'.join(result))

        # Run maf_join.py command using subprocess
        maf_join_output = f'{gene_name}_maf_join_output.maf'
        maf_join_cmd = f'python2.7 maf_join.py {maf_join_input} > {maf_join_output}'
        process = subprocess.run(maf_join_cmd, shell=True, check=True)

        if process.returncode == 0:
            print("maf_join.py command completed successfully and maf_join_output.maf is created")
        else:
            print(f"maf_join.py command failed with return code {process.returncode}")


def run_relax_analysis(maf_join_output, relax_input_fasta, relax_input_tree):
    fasta_lines = []
    species_names = []
    with open(maf_join_output, 'r') as maf_file:
        for line in maf_file:
            if line.startswith('a'):
                continue
            elif line.startswith('s'):
                fields = line.strip().split()
                species_name = fields[1].split('.')[0]  # Extract species name without the .chrX suffix
                sequence = fields[6]

                if sequence[-3:] in ['TAA', 'TGA', 'TAG', '---']:
                    sequence = sequence[:-3]

                fasta_line = f'>{species_name}\n{sequence}'
                fasta_lines.append(fasta_line)
                species_names.append(f'>{species_name}')

    fasta_content = '\n'.join(fasta_lines)

    # Write the FASTA content to a file
    with open(relax_input_fasta, 'w') as fasta_file:
        fasta_file.write(fasta_content)

    # Write the species names to a file
    species_names_str = ','.join(species_names).replace('>', '')

    # Prune using tree_doctor
    maf_join_cmd = f'tree_doctor -P {species_names_str} -n hg38.470way.nh > {relax_input_tree}'
    process = subprocess.run(maf_join_cmd, shell=True, check=True)

    if process.returncode == 0:
        print("Newick tree is created successfully")
    else:
        print(f"Could not create tree, failed with return code {process.returncode}")

    # Run Relax command
    relax_cmd = f'hyphy ./home/SharmaLab/hyphy/hyphy-analyses/RELAX-scan/RELAX.bf --alignment "{relax_input_fasta}" --tree "{relax_input_tree}" --reference hg38 > {gene_name}_hyphy.log  2>&1 & '
    process = subprocess.run(relax_cmd, shell=True, check=True)

    if process.returncode == 0:
        print("RELAX analysis completed successfully")
    else:
        print(f"RELAX analysis failed with return code {process.returncode}")     
 





 
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

            #maf_bb_file = run_maf_index(chr_number)
            run_maf_extract(gene_name, chr_number, exon_coords)
            
            maf_gene_seqs = {}
            gene_seqs = {} 
            query_seq = {}
            missing_bases_maf={}
            exon_strands ={}
            cesar_input_file = f'{gene_name}_cesar_input.fa'
            cesar_out_file= f'{gene_name}_cesar_output.fa'
            maf_join_output = f'{gene_name}_maf_join_output.maf'
            relax_input_tree = f'{gene_name}_relax_input.nh'
            relax_input_fasta=f'{gene_name}_relax_input.fa'
            for exon_number in exon_numbers:
                sequence_data = generate_coordinate_tables(gene_name)
                generate_sequence_dict(sequence_data, gene_name, exon_number)
            write_ms_file(cesar_input_file, query_seq)
            check_sequence_similarity(gene_seqs, maf_gene_seqs)
            compare_gene_lengths(gene_seqs, maf_gene_seqs, missing_bases_maf)
            write_fasta_files(gene_name, gene_strand, gene_seqs, cesar_input_file, exon_strands)
            create_MSA_file(cesar_input_file, gene_name)
            process_cesar_out_file(cesar_out_file, gene_name)
            run_relax_analysis(maf_join_output, relax_input_fasta, relax_input_tree)
            
            
           
            print('#############################################################################')
           

