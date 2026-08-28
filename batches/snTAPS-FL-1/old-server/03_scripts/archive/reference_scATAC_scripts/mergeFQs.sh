cat > mergeFQs.sh <<'EOF'
#!/bin/bash

INPUT_DIR=""
OUTPUT_DIR=""
THREADS=4

usage() {
    echo "用法: $0 -i <输入目录> -o <输出目录> [-t 线程数]"
    exit 1
}

while getopts "i:o:t:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        \?) usage ;;
    esac
done

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    usage
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "输入目录不存在: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "[$(date +'%Y-%m-%d %H:%M:%S')] 开始扫描输入目录: $INPUT_DIR"

find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" -type f \
    | awk -F'/' '{print $NF}' \
    | cut -d'.' -f1 \
    | sort -u > "${OUTPUT_DIR}/unique_samples.tmp"

SAMPLE_COUNT=$(wc -l < "${OUTPUT_DIR}/unique_samples.tmp")
echo "[$(date +'%Y-%m-%d %H:%M:%S')] 共发现 $SAMPLE_COUNT 个样本，开始并行合并 (线程数: $THREADS)"

cat "${OUTPUT_DIR}/unique_samples.tmp" | xargs -I {} -P "$THREADS" bash -c '
    sample="$1"
    in_dir="$2"
    out_dir="$3"

    echo "  -> 开始合并样本: ${sample}"

    find "${in_dir}" -maxdepth 1 -name "${sample}.*R1.fastq.gz" -type f -print0 \
        | xargs -0 cat > "${out_dir}/${sample}_R1.fastq.gz"

    find "${in_dir}" -maxdepth 1 -name "${sample}.*R2.fastq.gz" -type f -print0 \
        | xargs -0 cat > "${out_dir}/${sample}_R2.fastq.gz"
' _ "{}" "$INPUT_DIR" "$OUTPUT_DIR"

rm -f "${OUTPUT_DIR}/unique_samples.tmp"
echo "[$(date +'%Y-%m-%d %H:%M:%S')] merge 完成: $OUTPUT_DIR"
EOF

chmod +x mergeFQs.sh