Scripts:
assembly_stats.py - calculates summary statistics for an assembled genome. Input = a fasta file.
1snp.py - randomly subsamples 1 snp per locus from a vcf file. Input = filtered vcf file with ID column similar to ipyrad output (e.g. 'locus3_position24')
count_n_sequences.py - count the number of sequences that only contain Ns. Input = a fasta file
remove_nscaffs.py - remove scaffolds with only Ns. Input = a fasta file
structure_job_creator2.py - creates parameter and PBS files for running STRUCTURE on a cluster
summarize_vcf.py - calculates summary statistics from a vcf file - Input = a vcf file from ipyrad