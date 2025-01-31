#!/usr/bin/bash

# ----------------------------------------------------------------------------------------------------------------------
# Genome file preparation
cd soybean_data/genome_file
# merge soybean genome, mt, chlp and rhizobia genome to a single fatsa file
cat Gmax_275_v2.0.fa mtchlp.fa USDA110.fa > Gmax_USDA_org.fa

# build genome index
hisat2-build Gmax_USDA_org.fa Gmax_USDA110 -p 20 
hisat2-build Gmax_rRNA.fa Gmax_rRNA -p 20

# in order to use the R script to annotate the TSS on the genome, we need to convert the gff3 format using gffread.
gffread Gmax_275_Wm82.a2.v1.gene.gff3 -o Gmax_v2.0_changed.gff3
cd ../..
# ----------------------------------------------------------------------------------------------------------------------

# Preprocess the original sequencing data
cd soybean_data/stripe-seq_data
# the Script including the following steps using nodule data as example.

#remove the adapters using cutadapt
cutadapt \
    -a AGATCGGAAGAGCACACG \
    -a "G{100}" -n 2 -j 20 --minimum-length 50 \
    -o Nod.1.rmadaptor3.fastq \
    LCS8044_WXT_Nodule-RT01_R1.fastq.gz 
cutadapt \
    -g GACGCTCTTCCGATCT \
    -j 20 --minimum-length 50 -n 4 \
    -o Nod.1.rm5.fastq \
    Nod.1.rmadaptor3.fastq 

# select the NNNNNNNNTATAGGG contained reads, the script is from https://github.com/zentnerlab/STRIPE-seq
perl ../../selectReadsByPattern.pl -m NNNNNNNNTATAGGG -o Nod.1.pattern -r Nod.1.rm5.fastq 

# change the fastq to fasta format
fastx_collapser -i  Nod.1.pattern.fq -o Nod.1.uniq.fa

# remove the TATAGGG adapter
cutadapt -g TATAGGG -j 20 --minimum-length 50 -n 4 -o Nod.1.clean.fq Nod.1.uniq.fa

# remove the reads from rRNA.
hisat2 \
    -t -p 20 \
    --max-intronlen 10000  \
    -f -x ../genome_file/Gmax_rRNA \
    -U Nod.clean.fa \
    -S nod.rRNA.sam \
    --un Nod.keep.fasta
# remove the reads to combinative genome.
hisat2 \
    -t -p 20 \
    --max-intronlen 10000  \
    -f -x ../genome_file/Gmax_USDA110 \
    -U Nod.keep.fasta \
    -S Nod.hisat.sam

# selected the high quality reads mapped to the genome rather than organelles and rhizobium. 
# we also ignore the reads mapped to scaffolds in the genome.
samtools view -@ 20 -h -q 30 Nod.hisat.sam | perl -nlae 'print if !(/JX|NC|scaff/i)' | samtools sort - -@ 20 -o Nod.soy.hisat.bam
samtools index Nod.soy.hisat.bam -@ 20

cd ../..
# ----------------------------------------------------------------------------------------------------------------------
