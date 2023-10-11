#Author: Keshav Hibare, PhD student, University of Limerick
#Contact:  Keshav.hibare@gmail.com

#!/bin/bash

# Set the path to the ucsc_tools directory
ucsc_tools_path=~/ucsc_tools

# Loop over all directories in the current directory
for species_dir in */
do
  # Remove the trailing slash from the directory name to get the species name
  species="${species_dir%/}"

  # Generate the chrom.sizes file for this species
  $ucsc_tools_path/twoBitInfo "$species/$species.2bit" stdout | sort -k2rn > "$species/chrom.sizes"
done
