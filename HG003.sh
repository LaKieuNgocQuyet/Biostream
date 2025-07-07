#! /bin/bash
#>>>>>Variables<<<<<<
REFERENCE="/home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
REFERENCE_RTG=/home/dev_quyet/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.sdf
READ1="/home/dev_quyet/HG003/SRR2962692_1.fastq.gz"
READ2="/home/dev_quyet/HG003/SRR2962692_2.fastq.gz"
LOG="/home/dev_quyet/HG003/runtime.log"
BENCHMARK=/home/dev_quyet/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
BENCHMARK_SNP=/home/dev_quyet/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.SNPs.vcf.gz
BENCHMARK_INDEL=/home/dev_quyet/benchmark/HG003_GRCh38_1_22_v4.2.1_benchmark.Indels.vcf.gz
EXOME_BED=/home/dev_quyet/GRCh38/Twist_Exome_Core_Covered_Targets_hg38.bed

BWA_SAM=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.sam
BWA_BAM=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.bam
BWA_BAM_SORT=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.sort.bam
BWA_BAM_SORT_MARKED=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.sort.marked.bam
BWA_BAM_MATRIC=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.output.metrics.txt
BWA_BAM_SORT_MARKED_RECAL=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.sort.marked.recal.bam
BWA_BAM_BQSR_TABLE=/home/dev_quyet/HG003/BWA/Map/HG003_exome.BWA.recal_data.table
BWA_RAW_VCF_BCFTOOLS=/home/dev_quyet/HG003/BWA/VCF/HG003_exome.BWA.bcftools.raw.vcf.gz
BWA_RAW_VCF_FREEBAYES=/home/dev_quyet/HG003/BWA/VCF/HG003_exome.BWA.freebayes.raw.vcf
BWA_RAW_VCF_GATK=/home/dev_quyet/HG003/BWA/VCF/HG003_exome.BWA.gatk.raw.vcf.gz
BWA_RAW_VCF_BCFTOOLS_SNP=/home/dev_quyet/HG003/BWA/VCF/SNPs/HG003_exome.BWA.bcftools.raw_SNPs.vcf.gz
BWA_RAW_VCF_BCFTOOLS_INDEL=/home/dev_quyet/HG003/BWA/VCF/Indels/HG003_exome.BWA.bcftools.raw_Indels.vcf.gz
BWA_RAW_VCF_FREEBAYES_SNP=/home/dev_quyet/HG003/BWA/VCF/SNPs/HG003_exome.BWA.freebayes.raw_SNPs.vcf.gz
BWA_RAW_VCF_FREEBAYES_INDEL=/home/dev_quyet/HG003/BWA/VCF/Indels/HG003_exome.BWA.freebayes.raw_Indels.vcf.gz
BWA_RAW_VCF_GATK_SNP=/home/dev_quyet/HG003/BWA/VCF/SNPs/HG003_exome.BWA.gatk.raw_SNPs.vcf.gz
BWA_RAW_VCF_GATK_INDEL=/home/dev_quyet/HG003/BWA/VCF/Indels/HG003_exome.BWA.gatk.raw_Indels.vcf.gz

MINIMAP2_SAM=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.sam
MINIMAP2_BAM=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.bam
MINIMAP2_BAM_SORT=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.sort.bam
MINIMAP2_BAM_SORT_MARKED=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.sort.marked.bam
MINIMAP2_BAM_MATRIC=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.output.metrics.txt
MINIMAP2_BAM_SORT_MARKED_RECAL=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.sort.marked.recal.bam
MINIMAP2_BAM_BQSR_TABLE=/home/dev_quyet/HG003/Minimap2/Map/HG003_exome.Minimap2.recal_data.table
MINIMAP2_RAW_VCF_BCFTOOLS=/home/dev_quyet/HG003/Minimap2/VCF/HG003_exome.Minimap2.bcftools.raw.vcf.gz
MINIMAP2_RAW_VCF_FREEBAYES=/home/dev_quyet/HG003/Minimap2/VCF/HG003_exome.Minimap2.freebayes.raw.vcf
MINIMAP2_RAW_VCF_GATK=/home/dev_quyet/HG003/Minimap2/VCF/HG003_exome.Minimap2.gatk.raw.vcf.gz
MINIMAP2_RAW_VCF_BCFTOOLS_SNP=/home/dev_quyet/HG003/Minimap2/VCF/SNPs/HG003_exome.Minimap2.bcftools.raw_SNPs.vcf.gz
MINIMAP2_RAW_VCF_BCFTOOLS_INDEL=/home/dev_quyet/HG003/Minimap2/VCF/Indels/HG003_exome.Minimap2.bcftools.raw_Indels.vcf.gz
MINIMAP2_RAW_VCF_FREEBAYES_SNP=/home/dev_quyet/HG003/Minimap2/VCF/SNPs/HG003_exome.Minimap2.freebayes.raw_SNPs.vcf.gz
MINIMAP2_RAW_VCF_FREEBAYES_INDEL=/home/dev_quyet/HG003/Minimap2/VCF/Indels/HG003_exome.Minimap2.freebayes.raw_Indels.vcf.gz
MINIMAP2_RAW_VCF_GATK_SNP=/home/dev_quyet/HG003/Minimap2/VCF/SNPs/HG003_exome.Minimap2.gatk.raw_SNPs.vcf.gz
MINIMAP2_RAW_VCF_GATK_INDEL=/home/dev_quyet/HG003/Minimap2/VCF/Indels/HG003_exome.Minimap2.gatk.raw_Indels.vcf.gz
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/usr/bin/time -v \
fasterq-dump SRR2962692 --progress --split-files -O /home/dev_quyet/HG003 \
2>> "$LOG"
bgzip /home/dev_quyet/HG003/SRR2962692_1.fastq
bgzip /home/dev_quyet/HG003/SRR2962692_2.fastq

bwa index "${REFERENCE}"

/usr/bin/time -v \
    bwa mem -t 4 -R "@RG\tID:group1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1" \
    "${REFERENCE}" \
    "${READ1}" \
    "${READ2}" \
  > "${BWA_SAM}" \
2>> "$LOG"

/usr/bin/time -v \
samtools view -Sb "${BWA_SAM}" > "${BWA_BAM}" \
2>> "$LOG"

/usr/bin/time -v \
samtools sort -o "${BWA_BAM_SORT}" "${BWA_BAM}" \
2>> "$LOG"

/usr/bin/time -v \
gatk MarkDuplicates \
    -I "${BWA_BAM_SORT}" \
    -O "${BWA_BAM_SORT_MARKED}" \
    -M "${BWA_BAM_MATRIC}" \
2>> "$LOG"

/usr/bin/time -v bash -c "
  bcftools mpileup -Ou -f \"$REFERENCE\" \"$BWA_BAM_SORT_MARKED\" | \
  bcftools call -mv -Oz -o \"$BWA_RAW_VCF_BCFTOOLS\"
" 2>> "$LOG"

/usr/bin/time -v \
freebayes -f \
    "${REFERENCE}" \
    "${BWA_BAM_SORT_MARKED}" \
  > "${BWA_RAW_VCF_FREEBAYES}" \
2>> "$LOG"

bgzip "${BWA_RAW_VCF_FREEBAYES}"

/usr/bin/time -v \
gatk BaseRecalibrator \
  -I "${BWA_BAM_SORT_MARKED}" \
  -R "${REFERENCE}" \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_omni2.5.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/hapmap_3.3.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /home/dev_quyet/Known-sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O "${BWA_BAM_BQSR_TABLE}" \
2>> "$LOG"

/usr/bin/time -v \
gatk ApplyBQSR \
  -R "${REFERENCE}" \
  -I "${BWA_BAM_SORT_MARKED}" \
  --bqsr-recal-file "${BWA_BAM_BQSR_TABLE}" \
  -O "${BWA_BAM_SORT_MARKED_RECAL}" \
2>> "$LOG"

/usr/bin/time -v \
gatk HaplotypeCaller \
  -R "${REFERENCE}" \
  -I "${BWA_BAM_SORT_MARKED_RECAL}" \
  -O "${BWA_RAW_VCF_GATK}" \
2>> "$LOG"
##########################################################################################################
/usr/bin/time -v \
    minimap2 -x sr -a \
    "${REFERENCE}.mmi" \
    -R "@RG\tID:group1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:sample1" \
    "${READ1}" \
    "${READ2}" \
  > "${MINIMAP2_SAM}" \
2>> "$LOG"

/usr/bin/time -v \
samtools view -Sb "${MINIMAP2_SAM}" > "${MINIMAP2_BAM}" \
2>> "$LOG"

/usr/bin/time -v \
samtools sort -o "${MINIMAP2_BAM_SORT}" "${MINIMAP2_BAM}" \
2>> "$LOG"

/usr/bin/time -v \
gatk MarkDuplicates \
    -I "${MINIMAP2_BAM_SORT}" \
    -O "${MINIMAP2_BAM_SORT_MARKED}" \
    -M "${MINIMAP2_BAM_MATRIC}" \
2>> "$LOG"

/usr/bin/time -v \
gatk BaseRecalibrator \
  -I "${MINIMAP2_BAM_SORT_MARKED}" \
  -R "${REFERENCE}" \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_omni2.5.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/hapmap_3.3.hg38.vcf.gz \
  --known-sites /home/dev_quyet/Known-sites_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites /home/dev_quyet/Known-sites_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -O "${MINIMAP2_BAM_BQSR_TABLE}" \
2>> "$LOG"

/usr/bin/time -v \
gatk ApplyBQSR \
  -R "${REFERENCE}" \
  -I "${MINIMAP2_BAM_SORT_MARKED}" \
  --bqsr-recal-file "${MINIMAP2_BAM_BQSR_TABLE}" \
  -O "${MINIMAP2_BAM_SORT_MARKED_RECAL}" \
2>> "$LOG"

/usr/bin/time -v \
gatk HaplotypeCaller \
  -R "${REFERENCE}" \
  -I "${MINIMAP2_BAM_SORT_MARKED_RECAL}" \
  -O "${MINIMAP2_RAW_VCF_GATK}" \
2>> "$LOG"

/usr/bin/time -v \
freebayes -f \
    "${REFERENCE}" \
    "${MINIMAP2_BAM_SORT_MARKED}" \
  > "${MINIMAP2_RAW_VCF_FREEBAYES}" \
2>> "$LOG"

bgzip "${MINIMAP2_RAW_VCF_FREEBAYES}"

/usr/bin/time -v bash -c "
  bcftools mpileup -Ou -f \"$REFERENCE\" \"$MINIMAP2_BAM_SORT_MARKED\" | \
  bcftools call -mv -Oz -o \"$MINIMAP2_RAW_VCF_BCFTOOLS\"
" 2>> "$LOG"
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
tabix -p vcf "${BWA_RAW_VCF_BCFTOOLS}"
tabix -p vcf "${BWA_RAW_VCF_FREEBAYES}.gz"
tabix -p vcf "${MINIMAP2_RAW_VCF_BCFTOOLS}"
tabix -p vcf "${MINIMAP2_RAW_VCF_FREEBAYES}.gz"
####################################################################################################################
#BWA_bcftools
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_BCFTOOLS}" \
   --select-type-to-include SNP \
   -O "${BWA_RAW_VCF_BCFTOOLS_SNP}" 

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_BCFTOOLS}" \
   --select-type-to-include INDEL \
   -O "${BWA_RAW_VCF_BCFTOOLS_INDEL}"
#BWA_freebayes
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_FREEBAYES}.gz" \
   --select-type-to-include SNP \
   -O "${BWA_RAW_VCF_FREEBAYES_SNP}" 

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_FREEBAYES}.gz" \
   --select-type-to-include INDEL \
   -O "${BWA_RAW_VCF_FREEBAYES_INDEL}"
#BWA_gatk
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_GATK}" \
   --select-type-to-include SNP \
   -O "${BWA_RAW_VCF_GATK_SNP}"

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BWA_RAW_VCF_GATK}" \
   --select-type-to-include INDEL \
   -O "${BWA_RAW_VCF_GATK_INDEL}"
###############################################################################################################
#Minimap2_bcftools
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_BCFTOOLS}" \
   --select-type-to-include SNP \
   -O "${MINIMAP2_RAW_VCF_BCFTOOLS_SNP}"

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_BCFTOOLS}" \
   --select-type-to-include INDEL \
   -O "${MINIMAP2_RAW_VCF_BCFTOOLS_INDEL}"
#Minimap2_freebayes
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_FREEBAYES}.gz" \
   --select-type-to-include SNP \
   -O "${MINIMAP2_RAW_VCF_FREEBAYES_SNP}"

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_FREEBAYES}.gz" \
   --select-type-to-include INDEL \
   -O "${MINIMAP2_RAW_VCF_FREEBAYES_INDEL}"
#Minimap2_gatk
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_GATK}" \
   --select-type-to-include SNP \
   -O "${MINIMAP2_RAW_VCF_GATK_SNP}"

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${MINIMAP2_RAW_VCF_GATK}" \
   --select-type-to-include INDEL \
   -O "${MINIMAP2_RAW_VCF_GATK_INDEL}"
################################################################################################################
gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BENCHMARK}" \
   --select-type-to-include SNP \
   -O "${BENCHMARK_SNP}"

gatk SelectVariants \
   -R "${REFERENCE}" \
   -V "${BENCHMARK}" \
   --select-type-to-include INDEL \
   -O "${BENCHMARK_INDEL}"
###############################################################################################################
#RTG_SNPs
rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${BWA_RAW_VCF_BCFTOOLS_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/BWA_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${BWA_RAW_VCF_FREEBAYES_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/BWA_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${BWA_RAW_VCF_GATK_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/BWA_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${MINIMAP2_RAW_VCF_BCFTOOLS_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/Minimap2_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${MINIMAP2_RAW_VCF_FREEBAYES_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/Minimap2_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_SNP}" \
  -c "${MINIMAP2_RAW_VCF_GATK_SNP}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/SNPs/Minimap2_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=snp \
  -e "${EXOME_BED}"
#RTG_Indels
rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${BWA_RAW_VCF_BCFTOOLS_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/BWA_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${BWA_RAW_VCF_FREEBAYES_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/BWA_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${BWA_RAW_VCF_GATK_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/BWA_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${MINIMAP2_RAW_VCF_BCFTOOLS_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/Minimap2_bcftools \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${MINIMAP2_RAW_VCF_FREEBAYES_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/Minimap2_freebayes \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"

rtg vcfeval \
  -b "${BENCHMARK_INDEL}" \
  -c "${MINIMAP2_RAW_VCF_GATK_INDEL}" \
  -t "${REFERENCE_RTG}" \
  -o /home/dev_quyet/HG003/RTG/Indels/Minimap2_gatk \
  --ref-overlap \
  --all-records \
  --roc-subset=indel \
  -e "${EXOME_BED}"