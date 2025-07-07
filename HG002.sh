#! /bin/bash

wget -P /home/dev_quyet/GRCh38 https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz

wget -P /home/dev_quyet/GRCh38 https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/illumina-prep/exome/hg38_Twist_Bioscience_for_Illumina_Exome_2.5.bed
wget -P /home/dev_quyet/benchmark https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P /home/dev_quyet/benchmark https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget -P /home/dev_quyet/GRCh38 https://www.twistbioscience.com/sites/default/files/resources/2022-01/Twist_Exome_Core_Covered_Targets_hg38.bed

wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget -P /home/dev_quyet/Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

/usr/bin/time -v \
    bwa index /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
2> /home/dev_quyet/time.log

/usr/bin/time -v bash -c '
    bwa mem -t 4 -R "@RG\tID:group1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1" \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
    /home/dev_quyet/HG002/SRR2962669_1.fastq \
    /home/dev_quyet/HG002/SRR2962669_2.fastq | \
    samtools view -Sb - | \
    samtools sort -@ 4 -o /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.bam -' \
2>> /home/dev_quyet/time.log

gatk MarkDuplicates \
    -I /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.bam \
    -O /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.bam \
    -M /home/dev_quyet/HG002/BWA/Map/output.metrics.bwa.txt

/usr/bin/time -v bash -c \
'bcftools mpileup -Ou -f \
/home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
/home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.bam | \
bcftools call -mv -Oz -o /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf' \
2>> /home/dev_quyet/time.log

bgzip /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf
tabix -p vcf /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf.gz

gunzip /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
samtools faidx /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

/usr/bin/time -v \
freebayes -f \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.bam \
  > /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf \
2>> /home/dev_quyet/time.log

bgzip /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf
tabix -p vcf /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf.gz

gatk CreateSequenceDictionary \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -O /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict

gatk BaseRecalibrator \
  -I /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.bam \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_omni2.5.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/hapmap_3.3.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /home/dev_quyet/Known-sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O /home/dev_quyet/HG002/BWA/Map/HG002_exome.bwa.recal_data.table

gatk ApplyBQSR \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -I /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.bam \
  --bqsr-recal-file /home/dev_quyet/HG002/BWA/Map/HG002_exome.bwa.recal_data.table \
  -O /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.recal.bam

/usr/bin/time -v \
gatk HaplotypeCaller \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -I /home/dev_quyet/HG002/BWA/Map/HG002_exome.BWA.sort.marked.recal.bam \
  -O /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf \
2>> /home/dev_quyet/runtime.log

bgzip /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf
tabix -p vcf /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf.gz

rtg format -o /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta


rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/BWA/BWA_freebayes \
  --ref-overlap \
  --all-records \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/BWA/BWA_bcftools \
  --ref-overlap \
  --all-records \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/BWA/BWA_gatk \
  --ref-overlap \
  --all-records \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed
##########################################################################################################
/usr/bin/time -v \
    minimap2 -x sr -d \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.mmi \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
2>> /home/dev_quyet/time.log


/usr/bin/time -v \
    minimap2 -x sr -a \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.mmi \
    -R "@RG\tID:group1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1" \
    /home/dev_quyet/HG002/SRR2962669_1.fastq \
    /home/dev_quyet/HG002/SRR2962669_2.fastq \
  > /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sam \
2>> /home/dev_quyet/runtime.log

samtools view -Sb /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sam > /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.bam

samtools sort -o /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.bam /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.bam

gatk MarkDuplicates \
    -I /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.bam \
    -O /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.bam \
    -M /home/dev_quyet/HG002/Minimap2/Map/output.metrics.Minimap2.txt

gatk BaseRecalibrator \
  -I /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.bam \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_omni2.5.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/hapmap_3.3.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /home/dev_quyet/Known-sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.recal_data.table

gatk ApplyBQSR \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -I /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.bam \
  --bqsr-recal-file /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.recal_data.table \
  -O /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.recal.bam

/usr/bin/time -v \
gatk HaplotypeCaller \
  -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -I /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.recal.bam \
  -O /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.gatk.raw.vcf \
2>> /home/dev_quyet/time.log

/usr/bin/time -v \
freebayes -f \
    /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
    /home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.bam \
  > /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.freebayes.raw.vcf \
2>> /home/dev_quyet/time.log

/usr/bin/time -v bash -c \
'bcftools mpileup -Ou -f \
/home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
/home/dev_quyet/HG002/Minimap2/Map/HG002_exome.Minimap2.sort.marked.bam | \
bcftools call -mv -Oz -o /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.bcftools.raw.vcf' \
2>> /home/dev_quyet/time.log
####################################################################################################################
#BWA_bcftools
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf.gz \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.bcftools.raw_SNPs.vcf.gz 

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.bcftools.raw.vcf.gz \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.bcftools.raw_Indels.vcf.gz
#BWA_freebayes
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf.gz \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.freebayes.raw_SNPs.vcf.gz 

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.freebayes.raw.vcf.gz \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.freebayes.raw_Indels.vcf.gz
#BWA_gatk
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf.gz \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.gatk.raw_SNPs.vcf.gz 

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/BWA/VCF/HG002_exome.BWA.gatk.raw.vcf.gz \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.gatk.raw_Indels.vcf.gz
###############################################################################################################
#Minimap2_bcftools
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V HG002/Minimap2/VCF/HG002_exome.Minimap2.bcftools.raw.vcf \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.bcftools.raw_SNPs.vcf.gz

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.bcftools.raw.vcf \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.bcftools.raw_Indels.vcf.gz
#Minimap2_freebayes
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.freebayes.raw.vcf \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.freebayes.raw_SNPs.vcf.gz

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.freebayes.raw.vcf \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.freebayes.raw_Indels.vcf.gz
#Minimap2_gatk
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.gatk.raw.vcf \
   --select-type-to-include SNP \
   -O /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.gatk.raw_SNPs.vcf.gz

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/HG002/Minimap2/VCF/HG002_exome.Minimap2.gatk.raw.vcf \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.gatk.raw_Indels.vcf.gz
################################################################################################################
gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
   --select-type-to-include SNP \
   -O /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz

gatk SelectVariants \
   -R /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
   -V /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
   --select-type-to-include INDEL \
   -O /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz
###############################################################################################################
#RTG_SNPs
rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.bcftools.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/BWA_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.freebayes.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/BWA_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/SNPs/HG002_exome.BWA.gatk.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/BWA_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.bcftools.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/Minimap2_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.freebayes.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/Minimap2_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/SNPs/HG002_exome.Minimap2.gatk.raw_SNPs.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/SNPs/Minimap2_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed
#RTG_Indels
rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.bcftools.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/BWA_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.freebayes.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/BWA_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/BWA/VCF/Indels/HG002_exome.BWA.gatk.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/BWA_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.bcftools.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/Minimap2_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.freebayes.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/Minimap2_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

rtg vcfeval \
  -b /home/dev_quyet/benchmark/HG002_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz \
  -c /home/dev_quyet/HG002/Minimap2/VCF/Indels/HG002_exome.Minimap2.gatk.raw_Indels.vcf.gz \
  -t /home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf \
  -o /home/dev_quyet/HG002/RTG/Indels/Minimap2_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e /home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed