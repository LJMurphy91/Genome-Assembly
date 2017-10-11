#BACTERIAL READ ALIGNMENT, VARIANT CALLING, GENOME ASSMEMBLY AND ASSEMBLY QUALITY ASSESSMENT
#ALL FILES NEED TO BE IN THE SAME FOLDER 
#WORKS WITH PAIRED END READS

########################################################################################

#USE SAMTOOLS, BWA MEM, FREEBAYES AND SPADES

#NEED - REFERENCE FILE IN FASTA FORMAT, REFERENCE IN GFF FORMAT & SEQUENCING DATA IN FASTQ FORMAT POST QC ASSESSMENT

########################################################################################

#1. FORMAT REFERENCE SEQUENCE

samtools faidx <reference.fna>
bwa index <reference.fna>

for inputfile1 in <*1.fq>; do inputfile2=${inputfile1//1/2}; outfile=${i//fq/bam}; bwa mem -t <no. threads> - R '@RG\tID:<your ID>\tSM:<your ID>' <reference.fna> $inputfile1 $inputfile2 | samtools view -Shu - | samtools sort -o - x > $outfile; done

########################################################################################

#2. VARIANT CALLING

for bamfiles in *.bam; do outfile=${i//.bam/.vcf}; freebayes -f <reference.fna> --ploidy 1 $bamfiles > $outfile; done

########################################################################################

#3. BACKUP VCF FILE

for vcfFiles in *.vcfl do outfile=${vcfFiles//.vcf/_bkup.vcf}; cp $vcfFiles $outfiles; done

########################################################################################

#4. COMPRESS VCF FOR TABIX ACCESS

for vcfFiles in *.vcf; do bgzip $vcfFiles; done

########################################################################################

#5. TABIX FOR EASY ACCESS OF VCF FILES

for zippedvcf in *vcf.gz; do tabix -p zippedvcf; done

########################################################################################

#6. OVERLAY GENE FEATURES OVER SNPS TO IDENTIFY GENES WHICH CONTAIN SNPS

for zippedvcf in *bkup.vcf; do outfile=${zippedvcf//_bkup.vcf/_variants.txt}; bedtools intersect -wao -a <reference.gff> -b $i > $outfile; done

########################################################################################

#7. GENOME ASSEMBLY

for inputfile1 in <*1.fq>; do inputfile2=${inputfile1//1/2}; outdir=spades/; python spades.py -t <no. threads. --careful --only_assembler -1 $inputfile1 -2 $inputfile2 -o $outdir; done

########################################################################################

#8. QUAST ASSEMBLY ASSESSMENT

for assemblies in *spades/contigs.fasta; do python ~/quast.py -o / -R <reference.fna> -G <reference.gff> $assemblies; done




