#!/bin/bash 

set -euo pipefail

conda activate rnaseq

# --- Diretório onde estão os arquivos FASTQ originais
input_dir="$HOME/your/analysis"
concat_dir="$HOME/your/analysis/concat"
fastqc_dir="$HOME/your/analysis/fastqc"
fastp_dir="$HOME/your/analysis/fastp"

mkdir -p "$concat_dir"
mkdir -p "$fastp_dir"
mkdir -p "$fastqc_dir"


# --- Loop para encontrar e concatenar arquivos
for i in $(find "$input_dir" -type f -name '*_L00*.fastq.gz' | awk -F/ '{print $NF}' | awk -F'_L' '{print $1}' | sort -u); do
    cat "$input_dir"/"$i"_L00*_R1_*.fastq.gz > "$concat_dir"/"$i".R1.fastq.gz
    cat "$input_dir"/"$i"_L00*_R2_*.fastq.gz > "$concat_dir"/"$i".R2.fastq.gz
done

# --- FastQC em todos os arquivos concatenados
echo "Rodando FastQC..."
for file in "$concat_dir"/*.fastq.gz; do
    sample_name=$(basename "$file" .fastq.gz)

    # checa se já existe relatório
    if [[ -f "$fastqc_dir/${sample_name}_fastqc.html" || -f "$fastqc_dir/${sample_name}_fastqc.zip" ]]; then
        echo "FastQC já realizado para $sample_name, pulando..."
    else
        echo "----------------------------------------"
        echo "Processando: $sample_name"
        fastqc -o "$fastqc_dir" "$file"
    fi
done

# --- MultiQC para consolidar relatórios
echo "Rodando MultiQC..."
multiqc "$fastqc_dir" -o "$fastqc_dir"


# --- Filtar e Trimming para Paired-end
for file in "$concat_dir"/*.R*.fastq.gz; do
    sample_name=$(basename "$file" | sed 's/.R[12].fastq.gz//')
    
    echo "----------------------------------------"
    echo "Processando $sample_name..."

    fastp -f 1 -t 1 -3 -5 -W 4 -M 25 -l 40 \ #adapt for your dataset
        --in1 "$concat_dir"/"$sample_name".R1.fastq.gz \
        --in2 "$concat_dir"/"$sample_name".R2.fastq.gz \
        --out1 "$fastp_dir"/"$sample_name"_FP.R1.paired.fastq.gz \
        --out2 "$fastp_dir"/"$sample_name"_FP.R2.paired.fastq.gz \
        --unpaired1 "$fastp_dir"/"$sample_name"_FP.R1.unpaired.fastq.gz \
        --unpaired2 "$fastp_dir"/"$sample_name"_FP.R2.unpaired.fastq.gz \
        -h -R "$fastp_dir"/fastp.report

    echo "Finalizado para $sample_name."
done

# --- Filtros e Trimming para Single-end
for file in "$concat_dir"/*.fastq.gz; do
    sample_name=$(basename "$file" .fastq.gz)

    echo "----------------------------------------"
    echo "Processando $sample_name..."

    # se já tiver rodado fastp, pula
    if [[ -f "$fastp_dir/${sample_name}_FP.fastq.gz" ]]; then
        echo "fastp já rodado para $sample_name, pulando..."
        continue
    fi

    fastp -f 1 -t 1 -3 -5 -W 4 -M 25 -l 40 \
        --in1 "$file" \
        --out1 "$fastp_dir/${sample_name}_FP.fastq.gz" \
        -h "$fastp_dir/${sample_name}_fastp.html" \
        -j "$fastp_dir/${sample_name}_fastp.json"

    echo "Finalizado para $sample_name."
done

