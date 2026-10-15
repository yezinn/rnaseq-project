#!/usr/bin/env bash
# 전체 파이프라인 실행.  사용: bash run_all.sh
# - 단계별 시작/종료 시각: results/run_log.txt
# - 끝나거나 실패하면 macOS 알림센터로 알림
set -Eeuo pipefail
source config.sh
export RESULTS_DIR
mkdir -p "$RESULTS_DIR"
LOG="$RESULTS_DIR/run_log.txt"
# 터미널 출력 전체(에러 포함)를 파일에도 저장 → 터미널을 닫아도 에러 확인 가능
exec > >(tee -a "$RESULTS_DIR/full_output.log") 2>&1

notify() {  # macOS 알림 (다른 OS에서는 조용히 무시)
  command -v osascript >/dev/null 2>&1 && \
    osascript -e "display notification \"$2\" with title \"$1\" sound name \"Glass\"" || true
}
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

CURRENT="setup"
trap 'log "FAILED at step: $CURRENT"; notify "synecho-rnaseq 실패" "단계: $CURRENT — 터미널 확인"' ERR

if [[ ! -f samplesheet.csv ]]; then
  python scripts/00_fetch_metadata.py "$BIOPROJECT"
  echo
  echo ">>> samplesheet.csv 의 condition 열을 확인/수정한 뒤 다시 'bash run_all.sh' 실행하세요."
  exit 0
fi

step() {
  CURRENT="$1"; shift
  log "START $CURRENT"
  "$@"
  log "DONE  $CURRENT"
}

log "===== run start (SUBSAMPLE_READS=${SUBSAMPLE_READS}, STRAND=${STRAND}, RESULTS_DIR=${RESULTS_DIR}) ====="
step download   bash scripts/01_download.sh
step qc_trim    bash scripts/02_qc_trim.sh
step reference  bash scripts/03_reference.sh
step align      bash scripts/04_align.sh
step count      bash scripts/05_count.sh
# conda 환경의 R을 강제로 사용 (/usr/local/bin 등 시스템 R이 먼저 잡히는 문제 방지)
RSCRIPT="${CONDA_PREFIX:-}/bin/Rscript"
if [[ ! -x "$RSCRIPT" ]]; then
  CURRENT="deseq2 (conda 환경에 R 없음: $RSCRIPT)"
  false
fi
step deseq2     "$RSCRIPT" scripts/06_deseq2.R
log "===== ALL DONE ====="
notify "synecho-rnaseq 완료" "결과: ${RESULTS_DIR}/deseq2"
