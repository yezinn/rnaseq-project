# ===== 프로젝트 설정 (여기만 고치면 됨) =====

# GEO GSE103462 에 연결된 SRA BioProject
BIOPROJECT="PRJNA401742"

# 참조 게놈: Synechococcus elongatus PCC 7942 (RefSeq assembly)
# ⚠️ 실행 전 NCBI에서 accession이 맞는지 한 번 확인할 것
GENOME_ACC="GCF_000012525.1"

# 빠른 테스트용: 샘플당 사용할 리드 수 (0이면 전체 사용)
#   마감 전 첫 실행은 2000000(200만) 정도 권장 — Mac에서도 금방 끝남
SUBSAMPLE_READS="${SUBSAMPLE_READS:-2000000}"   # 실행 시 덮어쓰기 가능: SUBSAMPLE_READS=0 bash run_all.sh

# featureCounts 가닥 설정 (0=unstranded, 1=stranded, 2=reverse) — scripts/07_check_strand.sh 가 자동 갱신
STRAND=2

THREADS=4

# 디렉토리
DATA_DIR="data"
RAW_DIR="${DATA_DIR}/raw"
# 서브샘플 / 전체 리드 결과를 분리 저장 (서로 덮어쓰지 않도록)
if [[ "${SUBSAMPLE_READS}" -gt 0 ]]; then RUN_TAG=""; else RUN_TAG="_full"; fi
TRIM_DIR="${DATA_DIR}/trimmed${RUN_TAG}"
REF_DIR="${DATA_DIR}/reference"
RESULTS_DIR="results${RUN_TAG}"
QC_DIR="${RESULTS_DIR}/qc"
BAM_DIR="${RESULTS_DIR}/bam"
COUNT_DIR="${RESULTS_DIR}/counts"
DE_DIR="${RESULTS_DIR}/deseq2"
