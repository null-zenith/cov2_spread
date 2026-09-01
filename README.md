# cov2_spread
Scripts for processing big SARS-CoV-2 data used in paper "Global and country-level SARS-CoV-2 spread patterns inferred from near-identical genome pairs"

0. First, you need GISAID data (fasta file with SARS-CoV-2 genomes and tsv file with metadata). It is possible to do it with data from other sources, but you need to modify the code accordingly.  
1. Start with processing fasta file with genomes with Nextclade and obtaining the table output in tsv format.
2. Run cov2plus1_ptrns.R in command line while passing the path to Nextclade tsv table as nc_path. This script does the following: 1) filters short and low quality genomes; 2) splits all of the genomes with N substitutions to separate files called N.tsv (1.tsv, 2.tsv, etc.); 3) looks for pairs of genomes one substitution apart in the pairs of files N.tsv and N+1.tsv (1.tsv and 2.tsv, then 2.tsv and 3.tsv, etc.) and writes them into a file. You can skip to the part 3 if the previous steps were already completed previously, just in case you need it.
3. Next you need to prepare a bootstrap table that contains information about genomes and locations where they were observed. This table is tab-separated and have only two columns and each row represents a genome. In the first column ("substitutions") substitution pattern in kept in Nextclade format (comma-separated substitutions "C1234T,A2345G,T22222C"). To make these tables more lightweight, only the substitution patterns of "connected" genomes (identical genomes or a pair of genomes one substitution away from each other), other substitution patterns are changed to NA instead. The second column (name is not important, but let's call this column "subregion") is a territory where it was observed. The value of "subregion" columns can be a country, or dirtrict, or UN subregion, so you have to determine desirable level of detalization of geographical information (some countries have too few genomes and it is reasonable to merge them together, but if you want to explore within-country diversity of connections, you need to split it). Examples of this process are presented in cov2make_boot_tbls.R. Example of boot table is boot_tbl.example.tsv - a few first rows of the table used in the paper.
4. Run cov2bootstrap.R script. If boot_tbl is not passed as a parameter, connection is inferred ONLY from identical genomes.
5. Run cov2aggregate.R script.


# Examples

an example of command line for cov2plus1_ptrns.R
```
Rscript cov2plus1_ptrns.R --nc_path sequences.tsv
```


an example of command line for cov2bootstrap.R
```
cov2bootstrap.R --boot_tbl boot_tbl.example.tsv --p1m1_tbl p1m1_from_to.example.tsv --N 100 
```

# cov2INF_SIM

cov2INF_SIM.R is just a script that is run from the Rstudio, you need to pass the log-transformed RR table to it and run it manually to assess probability of importation inferred from RR table
