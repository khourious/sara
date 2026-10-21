#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------
# Pipeline para testar se o tamanho do arquivo é o esperado no download
# ----------------------------------------

# Caso tenha feito download de amostras fastq 
# use este pipeline para checar se o tamanoh esperado foi baixado efetivamente

# Para isso, recomendo que tenha feito um script usando o wget para downlaod (https://sra-explorer.info/)
# e aplicado um nohup para execução, pois assim gera um arquivo de log

# algo como
# nohup bash seu_script.sh &

# Variáveis
OUTDIR="path/to/your/download/files"
LOGFILE="$OUTDIR/check_sizes.log"
OUTTXT="path/to/your/log/output"

: > "$LOGFILE"

echo "Pipeline outubro/2026" | tee -a "$LOGFILE"
echo "_____________________" | tee -a "$LOGFILE"

echo "Checando tamanhos com base em $OUTTXT..." | tee -a "$LOGFILE"

# extrair pares arquivo|tamanho esperado do out.txt
grep "SIZE" "$OUTTXT" | while read -r line; do
  fname=$(echo "$line" | awk '{print $3}')
  expected=$(echo "$line" | awk '{print $5}')

  localpath="$OUTDIR/$fname"
  actual=0
  if [[ -f "$localpath" ]]; then
    actual=$(stat -c%s "$localpath")
  fi

  echo "----" | tee -a "$LOGFILE"
  echo "Arquivo: $fname" | tee -a "$LOGFILE"
  echo "Esperado: $expected bytes" | tee -a "$LOGFILE"
  echo "Local:    $actual bytes" | tee -a "$LOGFILE"

  if [[ "$expected" == "$actual" ]]; then
    echo "✅ OK: arquivo completo" | tee -a "$LOGFILE"
  else
    echo "⚠️ Problemas em $fname" | tee -a "$LOGFILE"
  fi
done

echo "Checagem concluída. Veja $LOGFILE para detalhes."
