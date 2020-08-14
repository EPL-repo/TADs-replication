

Prepare sam files using commands suitable for single-end sequencing, for example:

bwa mem -t28   allchr_masked.fasta    repl_file_R1.fastq  >   repl_file_R1.sam   (highest quality: 60)

or

bwa aln -t28   allchr_masked.fasta    repl_file_R1.fastq  >   repl_file_R1.sai
bwa samse -n 1 allchr_masked.fasta    repl_file_R1.sai  repl_file_R1.fastq  >  repl_file_R1.sam   (highest quality: 37)


There is no need for paired-end sequencing, when the bin size of the final analysis is 10 Kb.

The scripts look for codes obtained for single-end sequencing.

If you have paired-end data, process only the first read (R1), or,
otherwise, modify the script to recognize the paired-end codes.

VERY IMPORTANT: Data must be aligned to the human masked genome version: GRCh37/hg19.

-----

Use only the reads with the highest quality score (60 or 37) and with complete match to the genome. 


for EdUseq samples:

calculate_number_hits_per_bin_v1.pl  10000  37  100  repl_file_R1.sam

calculate_sigma_from_number_hits_per_bin_v1.pl  repl_file_R1_bin-size_10000_quality_37_chr1-X_adjbin_0b.csv

plot_sigma_values_v1.pl    gm  chr7  701  500  1 1     4.7 1

(for EdUseq samples, plot the smoothened data: variables 1 or 3)


for EUseq samples:

calculate_number_hits_per_bin_EUseq_v1.pl  10000  37  100  EUseq_file_R1.sam

calculate_sigma_from_number_hits_per_bin_EUseq_v1.pl   EUseq_file_R1--EUfor_bin-size_10000_quality_37_chr1-X_adjbin_EU_0b.csv
calculate_sigma_from_number_hits_per_bin_EUseq_v1.pl   EUseq_file_R1--EUrev_bin-size_10000_quality_37_chr1-X_adjbin_EU_0b.csv

plot_sigma_values_EUseq_v1.pl   lm  chr7  701  500    2 2 2 2     1 1 1 1

(for EUseq samples, plot the trimmed data: variables 2 or 4)

-----

General format of csv files:

25,000 lines correspond to 25,000 bins of 10 Kb each
Line 1: nts 1-9,999 = bin: 0
Line 2: nts 10,000-19,999 = bin: 1, etc
Each line contains data for 23 chromosomes (chr 1-22 plus chr X), in the format:
chr_1_variable_1, chr_1_variable_2, chr_1_variable_n,, chr_2_variable_1, chr_2_variable_2, chr_2_variable_n,, etc

-----

Variables in file: GeneD_OiRD_ReplT_interpolated__bin-size_10000_chr1-X_0b.csv

variable 1: replication timing: 
   eSC: early S replicating
   eSU: early S replicating, but some replication in mid S
   mSC: mid S replicating
   mSU: mid S replicating, but some replication in late S
   lS: late S replicating
   UND: replication time undefined
   REP: repetitive (masked) sequence

variable 2: type of replication domain:
   CNRD: genomic domain replicated from constitutive origin
   OiRD: genomic domain replicated from oncogene-induced origin
   UNKn: not known what type of origin this domain is replicated from (typically mid S or late S replication domains)

variable 3: gene annotation
   0: no gene in bin
   -50: one gene in bin
   higher negative number: more than one genes in bin

variable 4: gene annotation
   0: no intergenic sequences in bin
   -200: one intergenic element in bin
   higher negative number: more than one intergenic elements in bin (separated by one or more genes)

variable 5: gene annotation
   name of gene/intergenic elements
   examples:
     intergenic;n,
     DOC2B;r,   (DOC2B gene, has reverse orientation in the genome)
     FAM20C;f,  (FAM20 gene, has forward orientation in the genome)
     SH3YL1;r_intergenic;n_ACP1;f, (bin contains: SH3YL1 gene, intergenic element, ACP1 gene)
     ZNF595--ZNF718;g, (overlapping genes with different orientations)
     C9orf66--DOCK8;d, (overlapping genes with different orientations)

-----

Format of file: Origin_List_bin-size_10000_Table_S1.csv
   Columns and rows are the same as in Supplementary Information Table 1

 
 







