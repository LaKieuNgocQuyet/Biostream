FROM ubuntu
WORKDIR /home/seqflow
RUN mkdir /home/seqflow/sample && \
    mkdir /home/seqflow/reference && \
    mkdir /home/seqflow/known_sites && \
    mkdir /home/seqflow/output
COPY seqflow.sh /home/seqflow/seqflow.sh
RUN chmod +x /home/seqflow/seqflow.sh
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
RUN wget -P /home/Variant_caller https://github.com/broadinstitute/gatk/releases/download/4.6.2.0/gatk-4.6.2.0.zip && \
    unzip /home/Variant_caller/gatk-4.6.2.0.zip && \
    rm /home/Variant_caller/gatk-4.6.2.0.zip 
RUN python3 -m venv /home/Variant_caller/Myvenv
ENV PATH="/home/Variant_caller/gatk-4.6.2.0:$PATH"
ENV PATH="/home/Variant_caller/Myvenv/bin:$PATH"
