# cov2_spread
Scripts for processing big SARS-CoV-2 data used in paper "Global and country-level SARS-CoV-2 spread patterns inferred from near-identical genome pairs"

0. First, you need GISAID data (fasta file with SARS-CoV-2 genomes and tsv file with metadata). It is possible to do it with data from other sources, but you need to modify the code accordingly.  
1. Start with processing fasta file with genomes with Nextclade and obtaining the table output in tsv format.
2. Run cov2plus1_ptrns.R in command line while passing the path to Nextclade tsv table as nc_path. This script does the following: 1) filters short and low quality genomes; 2) splits all of the genomes with N substitutions to separate files called N.tsv (1.tsv, 2.tsv, etc.); 3) looks for pairs of genomes one substitution apart in the pairs of files N.tsv and N+1.tsv (1.tsv and 2.tsv, then 2.tsv and 3.tsv, etc.) and writes them into a file. You can skip to the part 3 if the previous steps were already completed previously, just in case you need it.
3. sdsdsd







# Examples

an example of command line for cov2plus1_ptrns.R
```
Rscript cov2plus1_ptrns.R --nc_path sequences.tsv
```
