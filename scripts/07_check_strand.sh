#!/usr/bin/env bash
# 라이브러리 가닥 특이성 판별: featureCounts -s 0/1/2 의 할당률을 비교해 STRAND를 정하고 config.sh에 기록.
#   예) unstranded면 s1≈s2≈s0/2,  stranded면 s1≈s0 이고 s2≈0,  reverse면 s2≈s0 이고 s1≈0
# 사용: bash scripts/07_check_strand.sh   (기존 서브샘플 BAM 사용 → 1분 이내)
set -euo pipefail
source scripts/_common.sh
out="$RESULTS_DIR/strand_check"
mkdir -p "$out"

bams=()
for run in $(runs); do bams+=("$BAM_DIR/${run}.bam"); done

echo -e "strand\tassigned\ttotal\trate" > "$out/strand_check.tsv"
for s in 0 1 2; do
  featureCounts -T "$THREADS" -s "$s" -F GTF -t gene -g locus_tag \
    -a "$REF_DIR/genome.gff" -o "$out/s${s}.txt" "${bams[@]}" > "$out/s${s}.log" 2>&1
  awk -v s="$s" 'NR>1{t=0; for(i=2;i<=NF;i++){tot[i]+=$i} if($1=="Assigned"){for(i=2;i<=NF;i++)a+=$i}}
       END{for(i in tot)T+=tot[i]; printf "%s\t%d\t%d\t%.4f\n", s, a, T, a/T}' \
       "$out/s${s}.txt.summary" >> "$out/strand_check.tsv"
done
column -t "$out/strand_check.tsv"

choice=$(awk 'NR>1{r[$1]=$4} END{
  if (r[1]/r[0] > 0.9) print 1; else if (r[2]/r[0] > 0.9) print 2; else print 0}' "$out/strand_check.tsv")
reason=$(awk -v c="$choice" 'NR>1{r[$1]=$4} END{printf "s1/s0=%.2f, s2/s0=%.2f -> STRAND=%s", r[1]/r[0], r[2]/r[0], c}' "$out/strand_check.tsv")
echo "$reason" | tee "$out/decision.txt"

sed -i.bak "s/^STRAND=.*/STRAND=${choice}/" config.sh && rm -f config.sh.bak
echo "[OK] config.sh 에 STRAND=${choice} 기록"
