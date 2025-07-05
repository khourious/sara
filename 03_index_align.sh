#!/bin/bash

set -euo pipefail

conda activate rnaseq

fastp_dir="$HOME/your/analisys/fastp"
genome_index="$HOME/your/analisys/hisat2_index/"
align_dir="$HOME/your/analisys/hisat2"
quant_dir="$HOME/your/analisys/stringtie"
rseq_qc="$HOME/your/analisys/rseqc/"

mkdir -p "$align_dir"
mkdir -p "$genome_index"
mkdir -p "$quant_dir"

# --- Download your genome and index
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.p14.genome.fa.gz -O "$genome_index"/GRCh38.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz -O "$genome_index"/gencode.v48.gtf.gz
unzip "$genome_index"/GRCh38.fa.gz
unzip "$genome_index"/gencode.v48.gtf.gz

# ---  Build HISAT2 index
hisat2-build "$genome_index"/GRCh38.fa "$genome_index"/GRCh38/ #The resulting index files will use 'GRCh38' as their prefix.

# --- Alinhar com HISAT2

for file in "$fastp_dir"/*_FP.R1.paired.fastq.gz; do
    # --- Extrair nome de amostras
    sample_name=$(basename "$file" | sed 's/_FP\.R1\.paired\.fastq\.gz//')
    sam_file="$align_dir/$sample_name.hisat2.sam"
    bam_file="$align_dir/$sample_name.hisat2.sorted.bam"

    # --- Verifica se já existe sam ou bam da amostra

    if [ -f "$sam_file" ] || [ -f "$bam_file" ]; then 
        echo "Alinhamento de $sample_name já existe. Pulando..."
        continue
    fi

    # --- Arquivos de entrada
    r1_file="$fastp_dir/${sample_name}_FP.R1.paired.fastq.gz"
    r2_file="$fastp_dir/${sample_name}_FP.R2.paired.fastq.gz"

    hisat2 --dta \
           -x "$genome_index/GRCh38" \
           -1 "$r1_file" \
           -2 "$r2_file" \
           -S "$sam_file" \
           --rna-strandness RF \
           --new-summary 

    echo "Alinhamento completado para $sample_name"
done

# --- Converter SAM a BAM

echo "Conversão de SAM a BAM"

for sam_file in "$align_dir"/*.hisat2.sam; do
    sample_name=$(basename "$sam_file" .hisat2.sam)
    bam_index="$bam_file.bai"
    
    # --- Elimina arquivos desnecessários
    if [ -f "$bam_file" ] && [ -f "$bam_index" ]; then
        echo "Archivo BAM $bam_file já existe. Eliminando SAM..."
        rm "$sam_file"
        continue
    fi

    echo "----------------------------------------"
    echo "Processando: $sample_name"

    # --- Converter
    echo "Convertendo SAM a BAM:"
    samtools view -b -o "${bam_file%.sorted.bam}.bam" "$sam_file"
    
    # --- Ordenar BAM 
    echo "Ordenando BAM:"
    samtools sort -o "$bam_file" "${bam_file%.sorted.bam}.bam"
    
    # --- Indexar BAM
    echo "Indexando BAM:"
    samtools index "$bam_file"

    # --- Validar o BAM gerado
    echo "Validando BAM:"
    samtools quickcheck "$bam_file" || { echo "Error: BAM corrompido em $sample_name"; exit 1; }

    # --- Eliminar arquivos temporários
    echo "Limpando arquivos temporários."
    rm "${bam_file%.sorted.bam}.bam"
    rm "$sam_file"

    echo "Conversão completada para $sample_name"
done

# --- Check strand specificity of the RNA-seq data
for bam_file in "$align_dir"/*.hisat2.sorted.bam; do
    sample_name=$(basename "$bam_file" .hisat2.sorted.bam) 
    infer_experiment.py -r "$genome_index"/gencode.v48.gtf -i "$bam_file" > "$rseq_qc"/"${sample_name}_strandness.txt"
done