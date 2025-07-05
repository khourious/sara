#!/bin/bash 

set -euo pipefail

conda activate rnaseq

# --- Diretório onde estão os arquivos FASTQ originais
input_dir="$HOME/your/analysis"
concat_dir="$HOME/your/analysis/concat"
fastp_dir="$HOME/your/analysis/fastp"

mkdir -p "$concat_dir"
mkdir -p "$fastp_dir"


# --- Loop para encontrar e concatenar arquivos
for i in $(find "$input_dir" -type f -name '*_L00*.fastq.gz' | awk -F/ '{print $NF}' | awk -F'_L' '{print $1}' | sort -u); do
    cat "$input_dir"/"$i"_L00*_R1_001.fastq.gz > "$concat_dir"/"$i".R1.fastq.gz
    cat "$input_dir"/"$i"_L00*_R2_001.fastq.gz > "$concat_dir"/"$i".R2.fastq.gz
done


# --- Filtar e Trimming
for file in "$input_dir"/*.R*.fastq.gz; do
    sample_name=$(basename "$file" | sed 's/.R[12].fastq.gz//')

    echo "Processando $sample_name..."

    fastp -f 1 -t 1 -3 -5 -W 4 -M 25 -l 40 \ #adapt for your dataset
        --in1 "$input_dir"/"$sample_name".R1.fastq.gz \
        --in2 "$input_dir"/"$sample_name".R2.fastq.gz \
        --out1 "$fastp_dir"/"$sample_name"_FP.R1.paired.fastq.gz \
        --out2 "$fastp_dir"/"$sample_name"_FP.R2.paired.fastq.gz \
        --unpaired1 "$fastp_dir"/"$sample_name"_FP.R1.unpaired.fastq.gz \
        --unpaired2 "$fastp_dir"/"$sample_name"_FP.R2.unpaired.fastq.gz \
        -h -R "$fastp_dir"/fastp.report

    echo "Finalizado para $sample_name."
done
