import re

def replace_refseq_with_uscs(fasta_file, mapping_file):
    # Read the mapping file and create a dictionary
    mapping = {}
    with open(mapping_file, 'r') as map_file:
        for line in map_file:
            if line.startswith("refseq"):
                continue
            columns = line.strip().split('\t')
            refseq = columns[0]
            uscs = columns[1]  # Replace the index from 3 to 1
            mapping[refseq] = uscs

    # Read the FASTA file and replace refseq identifiers with uscs identifiers
    updated_lines = []
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                match = re.match(r'>(\S+)', line)
                if match:
                    refseq_id = match.group(1)
                    uscs_id = mapping.get(refseq_id)
                    if uscs_id:
                        line = '>' + uscs_id + '\n'
            updated_lines.append(line)
    
    # Write the updated lines back to the FASTA file
    with open(fasta_file_out, 'w') as fasta:
        fasta.writelines(updated_lines)

# Specify the filenames
fasta_file = 'GCF_014898055.1.fa'
fasta_file_out = 'HLtalOcc1.fa'
mapping_file = 'GCF_014898055.1.chromAlias_modified.txt'

# Call the function to replace refseq with uscs identifiers
replace_refseq_with_uscs(fasta_file, mapping_file)
