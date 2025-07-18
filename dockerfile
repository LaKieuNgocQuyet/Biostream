FROM ubuntu
WORKDIR /home/Biostream
RUN mkdir /home/Biostream/sample && \
    mkdir /home/Biostream/reference && \
    mkdir /home/Biostream/known_sites && \
    mkdir /home/Biostream/output
COPY seqflow.sh /home/Biostream/seqflow.sh
RUN chmod +x /home/Biostream/seqflow.sh
RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-venv
RUN apt-get install -y \
    wget \
    unzip
RUN apt-get install -y \
    fastqc \
    trimmomatic \
    minimap2 \
    bwa \
    samtools \
    bcftools

RUN wget -P /home/Biostream https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip && \
    unzip /home/Biostream/gatk-4.6.2.0.zip && \
    rm /home/Variant_caller/gatk-4.6.2.0.zip 
RUN python3 -m venv /home/Biostream/Myvenv
ENV PATH="/home/Biostream/gatk-4.6.2.0:$PATH"
ENV PATH="/home/Biostream/Myvenv/bin:$PATH"
