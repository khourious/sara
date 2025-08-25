#!/bin/bash 

set -euo pipefail

conda activate rnaseq

# --- Diretório onde estão os arquivos FASTQ originais
input_dir="$HOME/your/analysis"
concat_dir="$HOME/your/analysis/concat"
fastqc_dir="$HOME/your/analysis/fastqc"
fastqc_lane_dir="$HOME/your/analysis/fastqc_lane"
multiqc_dir="$HOME/your/analysis/multiqc"
fastp_dir="$HOME/your/analysis/fastp"

mkdir -p "$concat_dir"
mkdir -p "$fastp_dir"
mkdir -p "$fastqc_dir"
mkdir -p "$fastqc_lane_dir"
mkdir -p "$multiqc_dir"

# --- Loop para rodar FastQC nos arquivos originais
for fq in $(find "$input_dir" -type f -name '*_L00*_R[12]_*.fastq.gz'); do
    fastqc "$fq" --outdir="$fastqc_lane_dir"
done


# --- Loop para encontrar e concatenar arquivos
for i in $(find "$input_dir" -type f -name '*_L00*.fastq.gz' | awk -F/ '{print $NF}' | awk -F'_L' '{print $1}' | sort -u); do
    cat "$input_dir"/"$i"_L00*_R1_*.fastq.gz > "$concat_dir"/"$i".R1.fastq.gz
    cat "$input_dir"/"$i"_L00*_R2_*.fastq.gz > "$concat_dir"/"$i".R2.fastq.gz
done

# --- Loop para rodar FastQC nos arquivos concatenados


for fq in "$concat_dir"/*.fastq.gz; do
    fastqc "$fq" --outdir="$fastqc_dir"
done

multiqc "$fastqc_dir" --filename "Relatorio_Final_MultiQC.html" -o "$multiqc_dir"

# --- Filtar e Trimming
for file in "$concat_dir"/*.R*.fastq.gz; do
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
