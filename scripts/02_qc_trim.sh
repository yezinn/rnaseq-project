#!/usr/bin/env bash
# FastQC(원본) → fastp 트리밍 → FastQC(트리밍 후) → MultiQC 요약
set -euo pipefail
source scripts/_common.sh
mkdir -p "$TRIM_DIR" "$QC_DIR/raw" "$QC_DIR/trimmed" "$QC_DIR/fastp"

# FastQC는 보조 QC(fastp가 자체 QC 리포트를 만듦) → 실패해도 파이프라인은 계속 진행
fastqc() {
  command fastqc "$@" || echo "[warn] FastQC 실패 (계속 진행): $*" | tee -a "$QC_DIR/fastqc_errors.log"
}

for run in $(runs); do
  read -ra in <<< "$(fastq_for_run "$run" "$RAW_DIR")"
  fastqc -q -t "$THREADS" --outdir "$QC_DIR/raw" "${in[@]}"

  if [[ ${#in[@]} -eq 2 ]]; then
    fastp -i "${in[0]}" -I "${in[1]}" \
          -o "$TRIM_DIR/${run}_1.fastq.gz" -O "$TRIM_DIR/${run}_2.fastq.gz" \
          --detect_adapter_for_pe -q 20 -l 25 -w "$THREADS" \
          -j "$QC_DIR/fastp/${run}.json" -h "$QC_DIR/fastp/${run}.html"
    fastqc -q -t "$THREADS" --outdir "$QC_DIR/trimmed" \
           "$TRIM_DIR/${run}_1.fastq.gz" "$TRIM_DIR/${run}_2.fastq.gz"
  else
    fastp -i "${in[0]}" -o "$TRIM_DIR/${run}.fastq.gz" \
          -q 20 -l 25 -w "$THREADS" \
          -j "$QC_DIR/fastp/${run}.json" -h "$QC_DIR/fastp/${run}.html"
    fastqc -q -t "$THREADS" --outdir "$QC_DIR/trimmed" "$TRIM_DIR/${run}.fastq.gz"
  fi
done

multiqc -q -f -o "$QC_DIR" "$QC_DIR" || echo "[warn] MultiQC 실패 (계속 진행)"
echo "[OK] QC report: $QC_DIR/multiqc_report.html"
