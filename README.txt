Scripts: 
1snp_by_id.py - randomly subsamples 1 snp per locus from a vcf file. Input = filtered vcf file with locus names in ID column
assembly_stats.py - calculates summary statistics for an assembled genome. Input = a fasta file
build_uce_reference.py - returns the longest sequence from a multifasta file with < 1000bp and <10% missing data. Input = a fasta file
contig_lengths.py - returns contig length statistics. Input = a genomic fasta file
count_n_sequences.py - count the number of sequences that only contain Ns. Input = a fasta file
get_alignment_data.py - get number of taxa, length of alignment, etc from a nested directory of alignment files. Input = name of output file
get_threshold.py - get GMYC threshold values from a nested directory of GMYC output files. Input = name of output file
gff2cds.py - creates a cds file from genome annotation data. Input = a gff3 formatted file, a functional annotation file produced by GOFEAT, a genome fasta file
produce_smc_files.py - produces smc++ input files. Input = a vcf file, the population name, sample names, and the path to the smc++ program
remove_chromosomes.py - removes SNPs by chromosome name. Input = a vcf file, a file listing the chromosome names to be removed, the name of the output file
remove_nscaffs.py - remove scaffolds with only Ns. Input = a fasta file
remove_scaffs.py - remove scaffolds from a fasta file. Input = a fasta file, a file listing which scaffolds to remove
split_scaffolds.py - split scaffolds based on GenBank contamination screens. Input = a fasta file and a file listing contamination coordinates
structure_job_creator2.py - creates parameter and PBS files for running STRUCTURE on a cluster
subsample_structure.py - subsample a specified number of SNPs from a STRUCTURE file. Input = a STRUCTURE-formatted file, the number of SNPs to sample, the output file name
subsample_vcf.py - subample a specified number of SNPs from a vcf file. Input = a vcf file, the number of SNPs to sample, the output file name
summarize_vcf.py - calculates summary statistics from a vcf file - Input = a vcf file from ipyrad