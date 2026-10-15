#!/usr/bin/env bash
# ENA가 느릴 때: NCBI SRA(prefetch)로 받아 FASTQ.gz로 변환 → data/raw 에 배치 + 무결성/리드 수 검증
# 사용: bash scripts/fetch_sra.sh SRR6004392 SRR6004398
set -euo pipefail
source config.sh
mkdir -p "$RAW_DIR" sra_cache

for run in "$@"; do
  echo "===== $run ====="
  # 프로젝트 폴더에서 이미 prefetch 한 경우(./SRRxxx/SRRxxx.sra) 그대로 재사용
  [[ -s "$run/$run.sra" && ! -e "sra_cache/$run" ]] && mv "$run" "sra_cache/$run"
  [[ -s "sra_cache/$run/$run.sra" ]] || prefetch -O sra_cache "$run"

  tmp="sra_cache/${run}_fastq"
  rm -rf "$tmp"; mkdir -p "$tmp"
  fasterq-dump --threads "$THREADS" --outdir "$tmp" "sra_cache/$run/$run.sra"
  fq="$tmp/$run.fastq"
  [[ -s "$fq" ]] || { echo "[ERROR] $fq 없음 (paired-end면 _1/_2 로 나옴)"; ls "$tmp"; exit 1; }

  # 리드 수 검증: ENA 메타데이터(runs.tsv)의 read_count 와 비교
  n=$(( $(wc -l < "$fq") / 4 ))
  expect=$(awk -F'\t' -v r="$run" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}$h["run_accession"]==r{print $h["read_count"]}' runs.tsv)
  echo "[check] 리드 수: $n (기대값 $expect)"
  [[ "$n" == "$expect" ]] || { echo "[ERROR] 리드 수 불일치"; exit 1; }

  echo "[gzip] 압축 중..."
  gzip -c "$fq" > "$RAW_DIR/$run.fastq.gz.tmp"
  gzip -t "$RAW_DIR/$run.fastq.gz.tmp"
  mv "$RAW_DIR/$run.fastq.gz.tmp" "$RAW_DIR/$run.fastq.gz"
  touch "$RAW_DIR/$run.fastq.gz.verified"
  rm -f "$RAW_DIR/$run.sub.fastq.gz"     # 손상 파일에서 만든 옛 서브샘플 폐기
  rm -rf "$tmp"
  echo "[OK] $RAW_DIR/$run.fastq.gz (verified)"
done
