fna_file = 'GCF_002288905.1_ASM228890v2_genomic.fna'

# Open the *.fna file and manipulate the first word of lines that start with '>'
with open(fna_file, 'r') as f_in, open(fna_file + ".modified", 'w') as f_out:
    for line in f_in:
        if line.startswith('>'):
            parts = line.split(' ', 1)
            first_word = parts[0][1:] # Remove the '>' character
            rest_of_line = parts[1]

            # Replace '.' with 'v' in the first word
            first_word_modified = first_word.replace('.', 'v')
            # Add 'chrUn_' prefix to the first word
            first_word_modified = 'chrUn_' + first_word_modified
            f_out.write(f">{first_word_modified} {rest_of_line}")
        else:
            f_out.write(line)
