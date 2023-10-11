# CASVIT-D
Comapartive analysis studies of Vitamin D

<--- Multiz470way maf files ---->

maf_download.sh

it downloads the ~2TB maf files of Multiz470 from the UCSC website and extarcts it and generate MD5 of the downloaded maf files and compare with the MD5 present in the website.

<---Check the quality of the Multiz470way--->

We gathered database of highly conseverd genes (6301 genes) that are expected to be present in any complete genome from BUSCO obtained their locus.

quality_cord.py

This makes a list of chromosome name, starnd, start and end cordinate for exon for every gene from the mutiz470 maf files.

aggregate_plot.py, plotter.py

aggregate_plot.py creates a csv file in such a manner tht if a gene has some erros it will get 0 if it doesnt it will get 1. From this CSV file we generated 2 different aggreate plots one with aggregate_plot.py and another with plotter.py for better respresenataiom

<---- genome download ----->

We curated the accession numbers of most of the genome assemblies that are used in mutiz470way manually and made a list.

genome_download.sh

This takes Species_accession_list.txt and downloads assemblies from NCBI website.

Genome_extract.sh

This extarcts the genomes from compressed file and coverts genomes from fa format to 2bit format (fatoTwoBit). creates its info file (2BitInfo).

NOTE: A handful of genomes were dowloaded manually since it was a tricky task to find all the genomes at one place.

Note: Some genome assembly had diffrent chromosome/scafold names than the names present in mutiz470way. ex: in the assembly it was represented byrefseqAccession/genbankAccession in the maf file it was reprsented in ucscStyleName.

check_chrom_names.py

To resolve the above issue. first we checked how many assmblies have different representaion than that of multiz470ways using check_chrom_names.py

genomeheader_edit_*.py

Once we got the list of assmblies and saw the representaion, we needed to create mutiple scrits do to solve this issue. all those scripts are in called genomeheader_edit_*.py

<---MASTER SCRIPT--->

COORD2RELAX.py

You can call this a masterscript. It takes cordindate file(hg38 gene/s) in the genepred format in the as input. It extarcts the maf file of the exons, on by one from the Multiz470way maf assembly according to cordinates (mafExtract feature). Form that exon maf it obtains the other species cordinates. Then it gets the sequesces of the each exons from their respective 2bit genome assmbies (2bitTofa). After that it creates paiwise alignment of each species with hg38 gene (CESAR). Then it coverts the pairwise alignment into mutiple sequence alignent. Then it looks for the abnormal/bad/uninformative sequences and removes them. Next, it extracts the names of the species remaining and removes the other species from the phylogentic tree (tree_doctor). then it runs HYPHY RELAX using mutiple sequecence aligmnet and tree as inputs.
