#Reference_sample
##HG002_Whole_exome_sequencing_Paired-End
fastq-dump SRR2962669 --split-files -O /home/dev_quyet/HG002
#Reference_genome
##GRCh38
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz
##hg19
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/hg19.p13.plusMT.no_alt_analysis_set.fa.gz
#Benchmark
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
#Known-sites_hg38
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
#Known-sites_hg19
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz

wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/dbsnp_138.hg19.vcf.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/dbsnp_138.hg19.vcf.idx.gz

wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz

wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz
wget http://ftp.cbi.pku.edu.cn/pub/mirror/GATK/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
