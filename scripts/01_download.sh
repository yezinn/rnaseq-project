#!/usr/bin/env bash
# runs.tsv 의 FASTQ를 ENA에서 내려받고, 설정 시 서브샘플링한다.
set -euo pipefail
source config.sh
mkdir -p "$RAW_DIR"

# 열 이름으로 run / fastq_ftp 를 뽑음 (빈 칸이 있어도 열이 밀리지 않게 awk 사용)
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
            {print $h["run_accession"] "\t" $h["fastq_ftp"]}' runs.tsv |
while IFS=$'\t' read -r run ftp; do
  # samplesheet.csv 에 남겨둔 run만 받음 (불필요한 샘플은 시트에서 지우면 됨)
  grep -q ",${run}," samplesheet.csv || continue
  IFS=';' read -ra urls <<< "$ftp"
  for u in "${urls[@]}"; do
    f="$RAW_DIR/$(basename "$u")"
    ok="$f.verified"   # 무결성 검사를 통과한 파일 표시 (다음 실행 때 검사 생략)
    for attempt in 1 2 3; do
      if [[ ! -s "$f" ]]; then
        echo "[download] $run → $f (시도 $attempt)"
        curl -fL --retry 3 -o "$f" "https://${u}"
        rm -f "$ok"
      fi
      [[ -f "$ok" ]] && break
      echo "[check] gzip 무결성 검사: $f"
      if gzip -t "$f" 2>/dev/null; then
        touch "$ok"; break
      fi
      echo "[warn] 손상된 파일 → 다시 받음: $f"
      rm -f "$f" "${f%.fastq.gz}.sub.fastq.gz"   # 손상 파일에서 만든 서브샘플도 폐기
    done
    [[ -f "$ok" ]] || { echo "[ERROR] 3회 시도 후에도 손상: $f"; exit 1; }
  done

  if [[ "${SUBSAMPLE_READS}" -gt 0 ]]; then
    for fq in "$RAW_DIR/${run}"*.fastq.gz; do
      [[ "$fq" == *".sub.fastq.gz" ]] && continue
      out="${fq%.fastq.gz}.sub.fastq.gz"
      # paired-end는 같은 seed(-s42)로 뽑아야 R1/R2 짝이 유지됨
      [[ -s "$out" ]] || seqtk sample -s42 "$fq" "$SUBSAMPLE_READS" | gzip > "$out"
    done
  fi
done
echo "[OK] download done"
