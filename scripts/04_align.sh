#!/usr/bin/env bash
# bowtie2 정렬 → 정렬·인덱싱 → 정렬률 요약
# 세균은 스플라이싱이 없으므로 HISAT2/STAR 대신 bowtie2로 충분하다.
set -euo pipefail
source scripts/_common.sh
mkdir -p "$BAM_DIR"
idx="$REF_DIR/genome"

for run in $(runs); do
  bam="$BAM_DIR/${run}.bam"
  read -ra in <<< "$(fastq_for_run "$run" "$TRIM_DIR")"
  # 입력(트리밍 FASTQ)보다 BAM이 새로울 때만 건너뜀 → 입력이 바뀌면 자동 재정렬
  [[ -s "$bam" && "$bam" -nt "${in[0]}" ]] && continue
  if [[ ${#in[@]} -eq 2 ]]; then
    reads=(-1 "${in[0]}" -2 "${in[1]}")
  else
    reads=(-U "${in[0]}")
  fi
  bowtie2 -p "$THREADS" -x "$idx" "${reads[@]}" 2> "$BAM_DIR/${run}.bowtie2.log" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
  samtools flagstat "$bam" > "$BAM_DIR/${run}.flagstat"
done

echo -e "run\toverall_alignment_rate" > "$BAM_DIR/alignment_summary.tsv"
for run in $(runs); do
  rate=$(grep "overall alignment rate" "$BAM_DIR/${run}.bowtie2.log" | awk '{print $1}')
  echo -e "${run}\t${rate}" >> "$BAM_DIR/alignment_summary.tsv"
done
cat "$BAM_DIR/alignment_summary.tsv"
