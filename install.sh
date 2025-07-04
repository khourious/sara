#!/bin/bash

set -euo pipefail

# Instala o Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Define diretórios
CONDA_DIR="$HOME/miniconda3"
CONDA_BIN="$CONDA_DIR/bin"

# Verifica se o diretório existe
if [ ! -d "$CONDA_BIN" ]; then
    echo "Diretório $CONDA_BIN não encontrado. Verifique a instalação do Miniconda."
    exit 1
fi

# Arquivos de configuração de shell
SHELL_RC="$HOME/.bashrc"
ZSH_RC="$HOME/.zshrc"

# Remove entradas antigas para evitar duplicação
sed -i '/miniconda3\/bin/d' "$SHELL_RC"
sed -i '/miniconda3\/bin/d' "$ZSH_RC"

# Adiciona export PATH ao final dos arquivos
echo "export PATH=\"$CONDA_BIN:\$PATH\"" >> "$SHELL_RC"
echo "export PATH=\"$CONDA_BIN:\$PATH\"" >> "$ZSH_RC"

# Aplica as mudanças se estiver usando bash ou zsh
if [[ "$SHELL" == */zsh ]]; then
    source "$ZSH_RC"
elif [[ "$SHELL" == */bash ]]; then
    source "$SHELL_RC"
fi

echo "PATH atualizado com $CONDA_BIN para $SHELL"


#install stringtie v2.2.1
wget https://github.com/gpertea/stringtie/releases/download/v2.2.1/stringtie-2.2.1.Linux_x86_64.tar.gz
tar -xvzf stringtie-2.2.1.Linux_x86_64.tar.gz
#stringtie-2.2.1.Linux_x86_64/stringtie

#create environment and install tools
conda create -n rnaseq python=3.9
conda activate rnaseq
conda install -c bioconda fastqc multiqc fastp hisat2 rseqc samtools stringtie

#check versions of installed tools
fastqc --version
multiqc --version
fastp --version
hisat2 --version
rseqc --version
samtools --version
stringtie --version

CORE="$HOME/path/to/your/analysis"

#create directories for data organization
mkdir -p "$CORE"/concat
mkdir -p "$CORE"/fastqc
mkdir -p "$CORE"/fastqc_lane
mkdir -p "$CORE"/multiqc
mkdir -p "$CORE"/fastp
mkdir -p "$CORE"/hisat2
mkdir -p "$CORE"/hisat2_index
mkdir -p "$CORE"/genome
mkdir -p "$CORE"/rseqc
mkdir -p "$CORE"/strigtie
mkdir -p "$CORE"/strigtie/output
