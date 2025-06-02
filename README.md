# SPARHA
Scripts for analyzing distances between 5'-ends of reads on the + and - strands

pair_reads_positive_distance.py

The script processes a bam file to find paired reads from the forward and reverse strands within a specified distance.
It extracts the positions of reads from both strands and then pairs them based on proximity using binary search. 
Each read from both strands is used only once.
The distance between reads is always calculated as a positive value and falls within the range of 0 to the specified maximum distance. 
The results, including paired reads and discarded reads are saved to excel files.
Additionally, the script calculates and writes statistics about the distances between paired reads.

pair_reads_plus_minus_distance.py

The script processes a bam file to identify all read pairs from the forward and reverse strands that are located within a specified maximum distance.
It extracts the positions of primary alignments separately for each strand and searches for nearby read pairs using binary search. 
Each read can participate in multiple pairs, allowing for all possible combinations within the defined distance. 
The distance is calculated as a signed value indicating direction (negative if the reverse strand read is downstream) but is constrained by the absolute maximum distance. 
The results, including read pair coordinates and distances, are saved to excel files. 
The script also computes and exports a distribution of distances between paired reads.
