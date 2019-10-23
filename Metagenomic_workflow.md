# Metagenomic_workflow_Guaymas

## Quality check using fastqc

`fastqc /PATH/TO/READ_1_1.fasta -o /OUTPUT_DIR`

- check results online


## Quality trimming

`module load bbmap-38.44/38.44`

`bbduk.sh in=/PATH/TO/READ1_1.fa in2=/PATH/TO/READ2_1.fa out=READ1_1.fastq out2=READ2_1.fastq qtrim=rl trimq=20 minlength=50`

 - reads were trimmed on both sides -> qtrim=rl
 - quality was chosen as Q20 -> trimq=20 
 - minimum length defines a threshold of which read length is removed -> minlength=50 (everything below removed)
 - repeat 1) and 2) for all reads


## Merge reads

Both runs to be analyzed were taken from the same library, indicated by the same barcode. Therefore those runs were merged for further analysis after quality trimming.

`cat PATH/TO/FASTQ/READ1_1.fastq PATH/TO/FASTQ/READ1_2.fastq > READ1_merged.fastq`

`cat PATH/TO/FASTQ/READ2_1.fastq PATH/TO/FASTQ/READ2_2.fastq > READ2_merged.fastq`


## Kmer spectrum

Filter reads you do not want to assemble by generating a k-mer spectrum with bbnorm 

`module load bbmap-38.44/38.44`

`bbnorm.sh -Xmx64g in=/PATH/TO/READ1_merged.fastq in2=/PATH/TO/READ2_merged.fastq  hist=READ_all_khist.txt peaks=READ_all_khist_peaks.txt threads=4 passes=1`

- plot computed kmer spectrum READ_all_khist.txt: x=depth, y=counts (logarithmic scale)
- showing kmer counts over readdepth gives you a hint on how often specific kmers were sequenced
- having a low depth, means kmers are probably rare, these are then removed to only keep sequences, which are of importance
- initial point of first peak defines kmer cutoff point, half it to not use any importand reads, e.g. initial point=50 -> highbindepth=25 (X)
- to compare use READ_all_khist_peaks.txt and compare raphical chosen initial peak point with computed peak.


## Kmer cutoff

Remove all reads with k-mer counts below a defined cutoff that you can define with highbindepth.

`bbnorm.sh -Xmx64g in=//PATH/TO/READ1_merged.fastq in2=/PATH/TO/READ2_merged.fastq highbindepth=X outhigh=READ1_merged_higherX.fastq outhigh2=READ2_merged_higherX.fastq passes=1 threads=4 overwrite=t`


## Assembly
Make a de novo assembly from merged reads with kmer cutoff.

`module load spades/3.9.0` 

`spades.py -o ./ --pe1-1 PATH/TO/READ1_merged_higherX.fastq --pe1-2 PATH/TO/READ2_merged_higherX.fastq --meta -t 8`

- using the flag meta, since using metagenomic reads
- 8 threads


## GC_coverage

Compute GC and coverage by mapping reads back to your de novo assembly using bbmap

`module load bbmap-38.44/38.44` 

`bbmap.sh ref=/PATH/TO/assembly/scaffolds.fasta nodisk in=/PATH/TO/READ1_merged_higherX.fastq in2=/PATH/TO/READ2_merged_higherX.fastq covstats=READ_coverage threads=4`


## rRNA via barrnap

- search for 16S rRNA genes within the assembly, using barrnap
- database needed is SILVA_SSU.noLSU.masked.trimmed.fasta
- use genome bin tools accessory script for a formatted file you can directly use within gbtools in R (from github gbtools)

`module load vsearch/1.9.7` 

`perl /PATH/TO/gbtools/accessory_scripts/get_ssu_for_genome_bin_tools.pl -d /PATH/TO/DATABASE/SILVA/132/SILVA_SSU.noLSU.masked.trimmed.fasta -c 4 -a /PATH/TO/assembly/scaffolds.fasta -o rRNA_READS`


## Use R to plot your GC-coverge plot

`setwd("/PATH/TO/DIRECTORY/")`

`install.packages("sp")`

`install.packages("plyr")`

`install.packages("devtools")`

`install.packages("/PATH/TO/gbtools_2.6.0.tar.gz",repos = NULL,type = "source")`

`library(gbtools)`


- load your data into R, here only coverage and the SSU were used, but further options are available, see gbtools.

`a <- gbt(covstats = "/PATH/TO/READ_coverage",
         ssu = "/PATH/TO/rRNA_READS.ssu.tab")`

- covstats: computed before with bbmap
- SSU: ssu accessory script gbtools

`plot(a)`

`plot(a,cutoff=2000,ssu=TRUE,textlabel=TRUE)`

`a.bin1 <- choosebin(a, slice = 1,num.points = 10,save=TRUE,file="/PATH/TO/a.bin1.contigs.list")`

`points(a.bin1,col = "red", slice=1)`

`a.bin2 <- choosebin(a,slice = 1,num.points = 10,save=TRUE,file="/PATH/TO/a.bin2.contigs.list")`

`points(a.bin2,col = "blue", slice=1)`

`a.bin3 <- choosebin(a,slice=1,num.points=10,save=TRUE,file="/PATH/TO/a.bin3.contigs.list")`

`points(a.bin3,col = "green", slice=1)`

- cutoff defines minimum contig length, the lower the length, the "cloudier" the plot might become
- careful in choosing a too high contig length, short contigs, with sequences you might be looking for could be removed.
- grey dots symbolize single contigs, contig length proportional to size
- Bins will be shown by clouds of contigs with the same Coverage and GC values -> Bins with highest coverage are most dominant in the assembly
- use "choosebin" to manually choose a bin from the assembly or load Bins from Metabat/Metawatt or other program in plot
- chosen bins are saved as "a.bin1.scaffold.list" -> gives you a list with all contigs, which are found in the bin
- you can choose as many bins as you need/want


## Extract contigs with sequences from assembly

A list of contigs was assigned to your bin, now you need to extract the sequences belonging to these contigs

`seqtk subseq /PATH/TO/assembly/scaffolds.fasta a.bin1.scaffolds.lst > Bin1.fasta`

- repeat for all Bins


## checkm for Bin statistics

`checkm lineage_wf /PATH/TO/BINS/ checkm/ -x fasta -t 4 --tab_table -f checkm_out`

- checkm qa for extended bin summary

`checkm qa -o 2 -f checkm_qa_out --tab_table -t 4 ./checkm/lineage.ms checkm/`

- creates a folder checkm with all output and a table in your current directory
- after check of Bin statistics for all bins, either refine bins (not needed here) OR compute your bin targeted reassembly


---------------------------------------------------------------------------------------------------------------
## Bin targeted reassembly

## Readmapping with Bins for Bin targeted reassembly

`module load bbmap-38.44/38.44`

`module load samtools-1.5/1.5` 
 

- here the bin is mapped back to the reads.
- I mapped my bin back to the qualitry trimmed reads without kmer cutoff, to not lose any sequences 

`bbmap.sh in=./PATH/TO/READ1_merged.fastq  in2=/PATH/TO/READ2_merged.fastq  minid=0.98 threads=4 outm=READ_98.sam ref=PATH/TO/EXTRACTED_BIN/Bin1.fasta`

- for bin targeted reassembly an identity of 98% was chosen for the first 3 readmappings, afterwards continued with 99%
- from readmapping one file (READ_98.sam) was computed, this now needs to be separated again in two fastq files, to continue with a bin targeted reassembly

`samtools view -S -b READ_98.sam > READ_98.bam`

`samtools bam2fq READ_98.bam > READ_98.fastq`

`cat READ_98.fastq | grep '^@.*/1$' -A 3 --no-group-separator > READ_98_1.fastq`

`cat READ_98.fastq | grep '^@.*/2$' -A 3 --no-group-separator > READ_98_2.fastq`


## Bin targetd reassembly

`module load spades/3.9.0` 

`spades.py -o./reassembly --pe1-1 /PATH/TO/READ_98_1.fastq --pe1-2 /PATH/TO/READ_98_1.fastq --meta -t 8`

a)Compute new GC coverage, as 6)

`module load bbmap-38.44/38.44`

`bbmap.sh ref=/PATH/TO/reassembly/scaffolds.fasta nodisk in=/PATH/TO/reads/READ_98_1.fastq  in2=/PATH/TO/reads/READ_98_2.fastq covstats=READ_98_coverage threads=4`

b)Compute new rRNA, as 7)

`module load vsearch/1.9.7` 

`perl /PATH/TO/gbtools/accessory_scripts/get_ssu_for_genome_bin_tools.pl -d /PATH/TO/DATABASE/SILVA/132/SILVA_SSU.noLSU.masked.trimmed.fasta -c 4 -a /PATH/TO/reassembly/scaffolds.fasta -o rRNA_READS_98`

- load data back into gbtools (R)
- covstats: computed before with bbmap
- SSU: ssu accessory script gbtools


`a <- gbt(covstats = "PATH/TO/READ_coverage",
         ssu = "PATH/TO/READ.ssu.tab")`

`plot(a)`

`plot(a,cutoff=2000,ssu=TRUE,textlabel=TRUE)`

- choose Bins

`plot(a)`

`plot(a,cutoff=2000,ssu=TRUE,textlabel=TRUE)`

`a.bin1 <- choosebin(a, slice = 1,num.points = 10,save=TRUE,file="/PATH/TO/a.bin1.contigs.list")`

`points(a.bin1,col = "red", slice=1)`

`a.bin2 <- choosebin(a,slice = 1,num.points = 10,save=TRUE,file="/PATH/TO/a.bin2.contigs.list")`

`points(a.bin2,col = "blue", slice=1)`

`a.bin3 <- choosebin(a,slice=1,num.points=10,save=TRUE,file="/PATH/TO/a.bin3.contigs.list")`

`points(a.bin3,col = "green", slice=1)`


## Extract contigs with sequences from assembly

- a list of contigs was assigned to your bin, now you need to extract the sequences belonging to these contigs

`seqtk subseq /PATH/TO/assembly/scaffolds.fasta a.bin1.scaffolds.lst > Bin1.fasta`


## checkm for Bin statistics

`checkm lineage_wf ./ checkm/ -x fasta -t 4 --tab_table -f checkm_out`

- checkm qa for extended bin summary

`checkm qa -o 2 -f checkm_qa_out --tab_table -t 4 ./checkm/lineage.ms checkm/`


## Annotation

- repeat from 12) until N50 value (checkM) does not get any higher. Use the final bin for annotation

### a) send to RAST

### b) PROKKA
- before using PROKKA headers have to be simplified, for this, using anvio

`module load anaconda3`

`source activate anvio5` 

`anvi-script-reformat-fasta` 

`module load Prokka/1.11`

`prokka --outdir ./Bin1/ --kingdom Archaea(or Bacteria) --gcode 11 --cpus 4 --rnammer /PATH/TO/BinX_short.fasta`  


--------------------------------------------------------------------------------------------------------------

## PhyloFlash 16S rRNA gene fragment abundance analysis

For first insights into the metagenomic reads

`phyloFlash.pl -dbhome /PATH/TO/DATABASE/SILVA_132/ -lib READ (name) -read1 /PATH/TO/TRIMMED/READ/READ1_1.fastq -read2 /PATH/TO/TRIMMED/READ/READ1_2.fastq -readlength 250 -id 90`

- readlength to set expected readlength
- id to set the minimum id (%), which is allowed

--------------------------------------------------------------------------------------------------------------

## Classification of bins

### a)GTDB_tk uses bacterial and archaeal marker genes to assign taxonomy to bins. 

`module load anaconda2/2.4.0` 

`module load prodigal/2.6.3` 

`module load hmmer/3.1b2` 

`module load pplacer18/1.1.alpha18` 

`module load FastANI/1.1.0` 

`module load FastTreeMP/2.8.1` 

`cd /PATH/TO/GTDB/OUTPUT/FOLDER`

`gtdbtk classify_wf --genome_dir PATH/TO/BINS/ --out_dir ./OUT --cpus 4 -x fasta`

- gives you a table with GTDB_tk classification, ANI (Average nucleotide identity) FastANI taxonomy, reference radius etc. 
- ANI is calculated with FastANI, which is included in the GTDB_tk classify workflow.


### b)BLAST 16 S rRNA genes found in bin in BLASTn suite https://blast.ncbi.nlm.nih.gov
- uses only 16S rRNA genes for taxonomy

## Find 16S rRNA genes within annotated Bin (Prokka)

- to extract 16S sequences from Prokka results use samtools
- look into .ffn file via less and search for 16S using /16S
- use header >PROKKA_XXXXX to extract single sequence

`module load samtools`

`samtools faidx PATH/TO/FILE.ffn PROKKA_XXXXX > FILE_16S.fasta`


## RNAmmer for 16S rRNA genes

Find 16 S rRNA genes within Bin, which is not annotated yet or Assembly-> RNAmmer/Barrnap (here RNAmmer)

`module load RNAmmer/1.2` 

`cd /PATH/TO/OUTPUT/FOLDER/RNAmmer`

`rnammer -S arc -m ssu -xml READ_all_arc.xml -gff READ_all_arc.gff  -h READ_all_arc.hmmreport -f READ_all_arc.fasta  < /PATH/TO/BIN/BinX.fasta``` 

`rnammer -S bac -m ssu -xml READ_all_bac.xml -gff READ_all_bac.gff  -h READ_all_bac.hmmreport -f READ_all_bac.fasta  < /PATH/TO/BIN/BinX.fasta`

- 16S fasta files can then be used for calculating a phylogenetic tree in ARB.
- same for other genes you are looking for, within .faa file (e.g. McrA)
--------------------------------------------------------------------------------------------------------------

## 16S Phylogeny Tree ARB

- loading sequences into database
- aligning sequences using SINA alignment and additionally manual alignment, if bases are not aligned correctly
- additionally more sequences were included in the database, e.g. 6215 more sequences of Euryarchaeota for the Bin tree
- add sequences to the already existing tree

- mark sequences from within the phylum/class, including neighbors from the sequences of the bin (about 200-300 in total)
- mark an outgroup
- calculate a filter over common aligned regions with minimum similarity of 50%, second filter: termini to restrict calculations to boundaries of 16S
- calculate maximum likelihood tree with RAxML 7 (DNA) algorithm and rapid bootstrap analysis 

- same for each Phylogeny tree in this Thesis
--------------------------------------------------------------------------------------------------------------

## Phylogenetic analysis of methyl-CoM reductase (mcrA)

#### a)SeaView to align the sequences to overall alignment
- use cat command to include sequences in overall alignment
- fit sequences in alignment
- exclude parts of sequences in the beginning and end that might be longer than other sequences
- Save as Fasta format ALIGNMENT.fasta

#### b)ZORRO to create a masking filter for the alignment

`zorro ALIGNMENT.fasta > YOUR_MASK.mask`

`awk '{printf "%.0f\n", $1}' YOUR_MASK.mask | cat > INTEGER_MASK.MASK`

`raxmlHPC-PTHREADS -m PROTGAMMALG -f a -N autoMRE -p RANDOM_NUMBER -s ALIGNMENT.fasta -a INTEGER_MASK.MASK  -k -x RANDOM_NUMBER -n OUTPUT -T THREADS`

- load file in iTOL (interactive Tree of Life) online and build McrA tree from it, by modifying branches
--------------------------------------------------------------------------------------------------------------
