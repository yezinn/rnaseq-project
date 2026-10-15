#!/usr/bin/env bash
# NCBI Datasets API로 참조 게놈(FASTA)과 주석(GFF)을 받고 bowtie2 인덱스를 만든다.
set -euo pipefail
source config.sh
mkdir -p "$REF_DIR"
cd "$REF_DIR"

if [[ ! -s genome.fna ]]; then
  curl -fL -o ref.zip \
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${GENOME_ACC}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
  unzip -o -q ref.zip -d ref_unzip
  cat ref_unzip/ncbi_dataset/data/"${GENOME_ACC}"/*.fna > genome.fna
  cp  ref_unzip/ncbi_dataset/data/"${GENOME_ACC}"/genomic.gff genome.gff
fi

echo "[info] 서열 목록 (염색체 + 플라스미드):"
grep '>' genome.fna

[[ -s genome.1.bt2 ]] || bowtie2-build -q --threads "$THREADS" genome.fna genome
samtools faidx genome.fna
echo "[OK] reference ready"
