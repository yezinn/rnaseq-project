#!/usr/bin/env bash
# featureCounts로 유전자별 리드 수를 센다 (locus_tag 기준, 예: Synpcc7942_0001).
set -euo pipefail
source scripts/_common.sh
mkdir -p "$COUNT_DIR"

bams=()
for run in $(runs); do bams+=("$BAM_DIR/${run}.bam"); done

# paired-end 여부: 첫 run 기준
pe_flag=()
read -ra first <<< "$(fastq_for_run "$(runs | head -1)" "$TRIM_DIR")"
[[ ${#first[@]} -eq 2 ]] && pe_flag=(-p --countReadPairs)

# 가닥 설정은 config.sh 의 STRAND (scripts/07_check_strand.sh 로 판별)
echo "[info] featureCounts -s ${STRAND}"
featureCounts -T "$THREADS" ${pe_flag[@]+"${pe_flag[@]}"} -s "$STRAND" \
  -F GTF -t gene -g locus_tag \
  -a "$REF_DIR/genome.gff" \
  -o "$COUNT_DIR/counts.txt" "${bams[@]}"

echo "[OK] $COUNT_DIR/counts.txt  (요약: counts.txt.summary)"
