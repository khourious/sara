
#!/bin/bash

set -euo pipefail

conda deactivate

fastp_dir="$HOME/your/analisys/fastp"
align_dir="$HOME/your/analisys/hisat2"
quant_dir="$HOME/your/analisys/stringtie"
prepde_script="$HOME/stringtie-2.2.1.Linux_x86_64/stringtie/prepDE.py3"

mkdir -p "$quant_dir"

# --- Quantificar com Stringtie

string_dir="$HOME"/stringtie-2.2.1.Linux_x86_64/stringtie

for bam_file in "$align_dir"/*.hisat2.sorted.bam; do
    sample_name=$(basename "$bam_file" .hisat2.sorted.bam)

    echo "Rodando StringTie v2.2.1 para $sample_name..."

    "$string_dir"/stringtie "$bam_file" \
              -G "$genome_index"/gencode.v48.gtf \
              -o "$quant_dir/$sample_name.gtf" \
              -A "$quant_dir/$sample_name.gene_abundance.tsv" \
              -e -B

    echo "StringTie v2.2.1 finalizado para $sample_name."
done

# --- Gerar Matriz de Contagem

# Gera o sample list
echo "Gerando sample_list.txt com separadores por TAB..."
> sample_list.csv
for gtf in "$quant_dir"/*.gtf; do
    sample=$(basename "$gtf" .gtf)
    echo -e "$sample\t$gtf" >> sample_list.csv
done

# Executa o prepDE.py
"$prepde_script" -i sample_list.csv

echo "Matriz de contagem gerada: gene_count_matrix.csv e transcript_count_matrix.csv"

mv gene_count_matrix.csv "$quant_dir"
mv gene_count_matrix.csv "$quant_dir"
mv sample_list.csv "$quant_dir"


#Rseqc strandness analysis
conda activate rnaseq

for bam_file in /home/joyce/hisat2/*.hisat2.sorted.bam; do
    sample_name=$(basename "$bam_file" .hisat2.sorted.bam) 
    infer_experiment.py -r /home/joyce/genome/gencode.v38.annotation.gtf -i "$bam_file" > /home/joyce/rseqc/"$sample_name"_strandness.txt
done