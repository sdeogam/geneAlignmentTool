# **Gene to Genome Aligner**

This program aligns a gene against a genome sequence using a heuristic algorithm inspired by basic alignment local search tool(BLAST). It is important to know the location of a gene on the genome sequence to get insights on its regulation. Also, knowledge of the location is very important in gene editing. 

The length of a genome is very large compared to the length of a gene (10^6 vs 10^3). This program takes the first 9 characters of a gene and searches them in the genome using a sliding window. The matches are the loci where the alignment begins.

The program uses a heuristic algorithm to begin aligning the gene against the genome, the heuristic being the length of the gene yet to be aligned. The scoring function is used which gives 0 to matches, 2 to mismatch and 3 to deletions. The search tree can also backtrack to look for the final state.

The program is written in python. Use this command to run.

'''python geneAlignmentTool.py'''

The program takes genome and gene sequences from fasta/fna files.

The genome file used in this example is [E. coli](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)
The gene used is [beta galactosidase](https://www.ncbi.nlm.nih.gov/gene/945006)
The input files can be changed by altering the path.

All the aligned sequences are stored in allAlignments.txt and the top 10 least scoring alignments are stored in top10Alignments.txt.
