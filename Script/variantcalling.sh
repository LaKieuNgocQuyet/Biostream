#! /bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Set_color <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
BLUE='\e[34m'
NC='\e[0m'
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Output_functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
function error () {
    echo -e "${RED}[ERROR]: $1${NC}"
    return
}
function info () {
    echo -e "${BLUE}[INFO]: $1${NC}"
    return 
}
function notice () {
    echo -e "${GREEN}[NOTICE]: $1${NC}"
    return
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
if [[ -z "$1" ]]; then
    error "input.sh file has been forgotten"
    echo "Example: bash seqflow.sh /home/seqflow/output/HG001_input.sh"
    exit 1
fi




INPUT="$1"
source "$INPUT"


info "=== Sample information ==="
info "SAMPLE_NAME = $SAMPLE_NAME"
info "READ_1 = $READ_1"
info "READ_2 = $READ_2"
info "REFERENCE_GENOME = $REFERENCE_GENOME"
info "OUTPUT_DIRECTORY = $OUTPUT_DIRECTORY"
info "KNOWN_SITES = ${KNOWN_SITES[*]}"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Directory_check <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if [[ ! -d "$OUTPUT_DIRECTORY" ]]; then
    echo "No OUTPUT_DIRECTORY could be found: $OUTPUT_DIRECTORY"
    exit 1
fi

# 3. Căn chỉnh reads
infor "Mapping/alignment"
minimap2 -x sr -a "$REFERENCE_GENOME".mmi "$READ_1" "$READ_2" \
    -R "@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:ILLUMINA" \
    > "$OUTPUT_DIRECTORY/$SAMPLE_NAME.aligned.sam"
samtools view -Sb "$OUTPUT_DIRECTORY/$SAMPLE_NAME.aligned.sam" > "$OUTPUT_DIRECTORY/$SAMPLE_NAME.aligned.bam"

# 4. Sort và index BAM
echo "Bước 2: Sort & Index BAM"
samtools sort -o "$OUTPUT_DIRECTORY/$SAMPLE_NAME.sorted.bam" "$OUTPUT_DIRECTORY/$SAMPLE_NAME.aligned.bam"
samtools index "$OUTPUT_DIRECTORY/$SAMPLE_NAME.sorted.bam"

# 5. (Tuỳ chọn) Hiệu chỉnh BQSR nếu có known sites
if [[ ${#KNOWN_SITES[@]} -gt 0 ]]; then
    echo "Bước 3: BQSR với known sites"
    gatk BaseRecalibrator \
        -R "$REFERENCE_GENOME" \
        -I "$OUTPUT_DIRECTORY/$SAMPLE_NAME.sorted.bam" \
        $(for ks in "${KNOWN_SITES[@]}"; do echo -n "--known-sites $ks "; done) \
        -O "$OUTPUT_DIRECTORY/$SAMPLE_NAME.recal.table"

    gatk ApplyBQSR \
        -R "$REFERENCE_GENOME" \
        -I "$OUTPUT_DIRECTORY/$SAMPLE_NAME.sorted.bam" \
        --bqsr-recal-file "$OUTPUT_DIRECTORY/$SAMPLE_NAME.recal.table" \
        -O "$OUTPUT_DIRECTORY/$SAMPLE_NAME.recal.bam"

    BAM_FINAL="$OUTPUT_DIRECTORY/$SAMPLE_NAME.recal.bam"
else
    echo "Không có known sites → bỏ qua BQSR"
    BAM_FINAL="$OUTPUT_DIRECTORY/$SAMPLE_NAME.sorted.bam"
fi

# 6. Gọi biến thể
echo "Bước 4: Gọi biến thể bằng GATK HaplotypeCaller"
gatk HaplotypeCaller \
    -R "$REFERENCE_GENOME" \
    -I "$BAM_FINAL" \
    -O "$OUTPUT_DIRECTORY/$SAMPLE_NAME.raw.vcf"

# 7. (Tuỳ chọn) Lọc biến thể nếu có ngưỡng lọc
echo "Bước 5: Lọc biến thể (nếu có)"
FILTER_OPTS=()
[[ ! -z "$FILTER_QUAL" ]] && FILTER_OPTS+=("--filter-expression" "QUAL < $FILTER_QUAL" --filter-name "LowQUAL")
[[ ! -z "$FILTER_QD" ]] && FILTER_OPTS+=("--filter-expression" "QD < $FILTER_QD" --filter-name "LowQD")
[[ ! -z "$FILTER_FS" ]] && FILTER_OPTS+=("--filter-expression" "FS > $FILTER_FS" --filter-name "HighFS")
[[ ! -z "$FILTER_SOR" ]] && FILTER_OPTS+=("--filter-expression" "SOR > $FILTER_SOR" --filter-name "HighSOR")
[[ ! -z "$FILTER_MQ" ]] && FILTER_OPTS+=("--filter-expression" "MQ < $FILTER_MQ" --filter-name "LowMQ")
[[ ! -z "$FILTER_MQRANKSUM" ]] && FILTER_OPTS+=("--filter-expression" "MQRankSum < $FILTER_MQRANKSUM" --filter-name "LowMQRankSum")
[[ ! -z "$FILTER_READPOSRANKSUM" ]] && FILTER_OPTS+=("--filter-expression" "ReadPosRankSum < $FILTER_READPOSRANKSUM" --filter-name "LowReadPosRankSum")

if [[ ${#FILTER_OPTS[@]} -gt 0 ]]; then
    gatk VariantFiltration \
        -R "$REFERENCE_GENOME" \
        -V "$OUTPUT_DIRECTORY/$SAMPLE_NAME.raw.vcf" \
        -O "$OUTPUT_DIRECTORY/$SAMPLE_NAME.filtered.vcf" \
        "${FILTER_OPTS[@]}"
    echo "Kết quả đã lọc: $SAMPLE_NAME.filtered.vcf"
else
    echo "Không có tiêu chí lọc. Giữ nguyên $SAMPLE_NAME.raw.vcf"
fi

echo "🎉 Pipeline hoàn tất cho mẫu: $SAMPLE_NAME"
