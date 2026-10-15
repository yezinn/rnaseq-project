#!/usr/bin/env python3
"""
ENA에서 BioProject의 run 목록을 받아 runs.tsv 와 samplesheet.csv 초안을 만든다.

- runs.tsv       : run별 FASTQ 다운로드 URL, single/paired 여부, 리드 수
- samplesheet.csv: sample, run, condition 열. condition은 sample_title에서
                   자동 추정하므로 반드시 눈으로 확인하고 고칠 것.

사용법: python scripts/00_fetch_metadata.py PRJNA401742
"""
import sys
import io
import requests
import pandas as pd

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ",".join([
    "run_accession", "sample_accession", "sample_title", "experiment_title",
    "library_layout", "library_strategy", "instrument_model",
    "read_count", "fastq_bytes", "fastq_ftp",
])


def fetch_runs(bioproject: str) -> pd.DataFrame:
    params = {"accession": bioproject, "result": "read_run",
              "fields": FIELDS, "format": "tsv"}
    r = requests.get(ENA_API, params=params, timeout=60)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    if df.empty:
        sys.exit(f"[ERROR] {bioproject}: run이 없습니다. accession을 확인하세요.")
    return df


def guess_condition(title: str) -> str:
    """sample_title 문자열로 조건을 추정 (GSE103462: WT vs rel-)."""
    t = str(title).lower()
    if "rel" in t and ("-" in t or "mut" in t or "delta" in t or "Δ" in t):
        return "rel_mutant"
    if "wt" in t or "wild" in t:
        return "WT"
    return "UNKNOWN"


def main():
    bioproject = sys.argv[1] if len(sys.argv) > 1 else "PRJNA401742"
    df = fetch_runs(bioproject)
    df = df[df["library_strategy"].astype(str).str.upper() == "RNA-SEQ"] \
        if "library_strategy" in df else df
    df.to_csv("runs.tsv", sep="\t", index=False)

    sheet = pd.DataFrame({
        "sample": df["sample_title"].astype(str)
                    .str.replace(r"[^A-Za-z0-9_]+", "_", regex=True).str.strip("_"),
        "run": df["run_accession"],
        "condition": df["sample_title"].map(guess_condition),
        "title": df["sample_title"],
    })
    sheet.to_csv("samplesheet.csv", index=False)

    print(f"[OK] run {len(df)}개 → runs.tsv, samplesheet.csv")
    print(sheet.to_string(index=False))
    print("\nlibrary_layout:", df["library_layout"].value_counts().to_dict())
    if (sheet["condition"] == "UNKNOWN").any():
        print("\n[!] condition이 UNKNOWN인 샘플이 있습니다. "
              "samplesheet.csv를 열어 직접 WT / rel_mutant 로 채우세요.")
    print("[!] 자동 추정 결과도 GEO 샘플 설명과 대조해 반드시 확인하세요.")


if __name__ == "__main__":
    main()
