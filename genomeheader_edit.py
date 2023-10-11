#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com

import json
import os
import sys
import glob

def extract_fields(dirs: str):
    output = "refseqAccession_genbankAccession_ucscStyleName.txt"

    for dir_ in dirs:
        jsonl_input = os.path.join(dir_, "sequence_report.jsonl")
        output_path = os.path.join(dir_, output)

        with open(jsonl_input, 'r') as f_in, open(output_path, 'w') as f_out:
            for line in f_in:
                data = json.loads(line)
                refseqAccession = data.get('refseqAccession')
                genbankAccession = data.get('genbankAccession')
                ucscStyleName = data.get('ucscStyleName')

                if refseqAccession is not None or genbankAccession is not None or ucscStyleName is not None:
                    if refseqAccession is not None:
                        refseqAccession_parts = refseqAccession.split('.')
                        f_out.write(f"{refseqAccession_parts[0]}\t")
                    else:
                        f_out.write("\t")

                    if genbankAccession is not None:
                        genbankAccession_parts = genbankAccession.split('.')
                        f_out.write(f"{genbankAccession_parts[0]}\t")
                    else:
                        f_out.write("\t")

                    if ucscStyleName is not None:
                        f_out.write(f"{ucscStyleName}\n")
                    else:
                        f_out.write("\n")

        # Call the manipulate_fna function for the current directory
        refseq_genbank_ucsc_file = output_path
        manipulate_fna(dir_, refseq_genbank_ucsc_file)

def manipulate_fna(dir_: str, refseq_genbank_ucsc_file: str):
    # Find the *.fna file in the specified directory
    fna_files = glob.glob(os.path.join(dir_, "*.fna"))
    if len(fna_files) == 0:
        print(f"No *.fna file found in directory {dir_}")
        return
    elif len(fna_files) > 1:
        print(f"Multiple *.fna files found in directory {dir_}")
        return
    fna_file = fna_files[0]

    # Read the refseqAccession_genbankAccession_ucscStyleName.txt file into a dictionary
    refseq_genbank_ucsc = {}
    with open(refseq_genbank_ucsc_file, 'r') as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) == 3:
                refseqAccession, genbankAccession, ucscStyleName = parts
                refseq_genbank_ucsc[refseqAccession] = ucscStyleName
                refseq_genbank_ucsc[genbankAccession] = ucscStyleName

    # Open the *.fna file and manipulate the first word of lines that start with '>'
    with open(fna_file, 'r') as f_in, open(fna_file + ".modified", 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                parts = line.split(' ', 1)
                first_word = parts[0][1:] # Remove the '>' character
                rest_of_line = parts[1]

                # Remove the version number from the first word
                first_word_parts = first_word.split('.')
                first_word_without_version = first_word_parts[0]

                if first_word_without_version in refseq_genbank_ucsc:
                    new_first_word = refseq_genbank_ucsc[first_word_without_version]
                    f_out.write(f">{new_first_word} {rest_of_line}")
                else:
                    f_out.write(f">{first_word_without_version} {rest_of_line}")
            else:
                f_out.write(line)



if __name__ == "__main__":
    dirs = sys.argv[1:]
    extract_fields(dirs)
