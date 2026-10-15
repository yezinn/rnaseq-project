#!/usr/bin/env bash
# 공통 함수: run의 입력 FASTQ 찾기 (서브샘플 여부 반영, SE/PE 자동 판별)
source config.sh

fastq_for_run() {
  # 사용: fastq_for_run <run> <dir> → R1 [R2] 경로 출력
  local run="$1" dir="$2" suffix=".fastq.gz"
  [[ "$dir" == "$RAW_DIR" && "${SUBSAMPLE_READS}" -gt 0 ]] && suffix=".sub.fastq.gz"
  if [[ -s "$dir/${run}_1${suffix}" ]]; then
    echo "$dir/${run}_1${suffix} $dir/${run}_2${suffix}"
  else
    echo "$dir/${run}${suffix}"
  fi
}

runs() { tail -n +2 samplesheet.csv | cut -d, -f2; }
