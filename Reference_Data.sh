#! /bin/bash
#Reference_sample
##HG002_Whole_exome_sequencing_Paired-End
mkdir -p ./HG002 
fastq-dump SRR2962669 --split-files -O ./HG002
#Reference_genome
##GRCh38
mkdir -p ./GRCh38
wget -P ./GRCh38 https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
##hg19
mkdir -p ./hg19
wget -P ./hg19 https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz
#Benchmark
mkdir -p ./Benchmark
wget -P ./Benchmark https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -P ./Benchmark https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

#==================================================================================================
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Known_sites_hg38 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#==================================================================================================
mkdir -p ./Known-sites_hg38
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget -P ./Known-sites_hg38 https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

#==================================================================================================
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Known_sites_hg19 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#==================================================================================================
mkdir -p ./Known-sites_hg19
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz

wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/dbsnp_138.hg19.vcf.gz
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/dbsnp_138.hg19.vcf.idx.gz

wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz

wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz

wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz
wget -P ./Known-sites_hg19 http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
