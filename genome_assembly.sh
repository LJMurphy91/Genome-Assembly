#BACTERIAL READ ALIGNMENT, VARIANT CALLING, GENOME ASSMEMBLY AND ASSEMBLY QUALITY ASSESSMENT
#ALL FILES NEED TO BE IN THE SAME FOLDER 
#WORKS WITH PAIRED END READS

########################################################################################

#USE KRAKEN, BBSPLIT, SAMTOOLS, BWA MEM, FREEBAYES AND SPADES

#NEED - REFERENCE FILE IN FASTA FORMAT, REFERENCE IN GFF FORMAT & SEQUENCING DATA IN FASTQ FORMAT POST QC ASSESSMENT

########################################################################################
#1. IDENTIFY CONTAMINATION IN SEQUENCING FILES AND REMOVE IT, USING KRAKEN AND BBSPLIT

for fastqFile1 in *1.fq; do fastqFile2=${fastqFile1//1/2}; outfile=${i//fq/out}; ~/kraken -db /path/to/DB --threads NUM --paired --fastq-input --output $outfile; done

#1A. CONVERT KRAKEN OUTPUT TO READABLE FILE
for krakenOut in *.out; do outfile=${krakenOut//out/report}; kraken-report --db /path/to/DB $krakenOut $outfile;done

#1B. DOWNLOAD FASTA FILES FROM GENBANK FOR EACH OF THE SPECIES OUTSIDE THE MYCOBACTERIUM CLADE WHICH COVER MORE THAN 0.01% OF READS

#1C. USE BBSPLIT TO REMOVE READS THAT MAP TO THESE "CONTAMINATION" GENOMES AND SAVE "CLEANED" READS TO NEW FILES

for fastqFile1 in *1.fq; do fastqFile2==${fastqFile1//1/2}; outFile1=${fastqFile1//.fq/\_cleaned.fq}; outFile2=${fastqFile2//.fq/\_cleaned.fq}; ~/bbsplit.sh in1=fastqFile1 in2=fastqFile2 ref=<comma separated list of species> basename=out_%.fq outu1=outFile1 outu2=outFile2

########################################################################################

#2. FORMAT REFERENCE SEQUENCE USING SAMTOOLS AND BWA INDEX, AND ALIGN READS TO REFERENCE USING BWA MEM

samtools faidx <reference.fna>
bwa index <reference.fna>

for inputfile1 in <*1_cleaned.fq>; do inputfile2=${inputfile1//1/2}; outfile=${i//fq/bam}; bwa mem -t <no. threads> - R '@RG\tID:<your ID>\tSM:<your ID>' <reference.fna> $inputfile1 $inputfile2 | samtools view -Shu - | samtools sort -o - x > $outfile; done

########################################################################################

#3. VARIANT CALLING USING FREEBAYES

for bamfiles in *.bam; do outfile=${i//.bam/.vcf}; freebayes -f <reference.fna> --ploidy 1 $bamfiles > $outfile; done

########################################################################################

#4. BACKUP VCF FILE

for vcfFiles in *.vcfl do outfile=${vcfFiles//.vcf/_bkup.vcf}; cp $vcfFiles $outfiles; done

########################################################################################

#5. COMPRESS VCF FOR TABIX ACCESS

for vcfFiles in *.vcf; do bgzip $vcfFiles; done

########################################################################################

#6. TABIX FOR EASY ACCESS OF VCF FILES

for zippedvcf in *vcf.gz; do tabix -p zippedvcf; done

########################################################################################

#7. OVERLAY GENE FEATURES OVER SNPS TO IDENTIFY GENES WHICH CONTAIN SNPS

for zippedvcf in *bkup.vcf; do outfile=${zippedvcf//_bkup.vcf/_variants.txt}; bedtools intersect -wao -a <reference.gff> -b $i > $outfile; done

########################################################################################

#8. GENOME ASSEMBLY USING SPADES

for inputfile1 in <*1.fq>; do inputfile2=${inputfile1//1/2}; outdir=${inputfile1//\_1.fq/}; python spades.py -t <no. threads. --careful --only_assembler -1 $inputfile1 -2 $inputfile2 -o $outdir; done

########################################################################################

#9. QUAST ASSEMBLY ASSESSMENT

for assemblies in *spades/contigs.fasta; do python ~/quast.py -o / -R <reference.fna> -G <reference.gff> $assemblies; done

########################################################################################

#10. GENOME ANNOTATION USING PROKKA

for folders in */; do cd $folders; prokka --locustag $i --genus GENUS --species SPECIES --strain $folders --gcode 11 --outdir prokka_$folders --usegenus contigs.fasta --force --cpus NUM; done



