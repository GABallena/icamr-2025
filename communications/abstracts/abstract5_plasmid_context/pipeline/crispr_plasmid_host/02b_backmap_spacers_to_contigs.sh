#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

REL="relabeled_contigs"
SPALL="work/spacers_fa/ALL.spacers.fasta"
OUTDIR="work/spacers_fa"
MAPALL="${OUTDIR}/spacer_to_contig.tsv"

command -v blastn >/dev/null || { echo "[ERR] blastn not found"; exit 1; }
[[ -s "$SPALL" ]] || { echo "[ERR] Missing $SPALL"; exit 1; }
mkdir -p "$OUTDIR"

# make per-sample spacer FASTAs
echo -n > "$MAPALL"  # truncate combined map
for fa in ${REL}/*_contigs.fa; do
  samp=$(basename "$fa" .fa)           # e.g., SAMPLE-01_S1_contigs
  sfa="${OUTDIR}/${samp}.spacers.fa"
  db="${fa}"                           # makeblastdb will add suffixes next
  tsv="${OUTDIR}/${samp}.spacer2contig.tsv"

  # Extract spacers for this sample (headers start with >SAMPLE|)
  awk -v s="${samp}|" '
    BEGIN{keep=0}
    /^>/ {keep = index($0, ">" s)==1}
    { if(keep) print }
  ' "$SPALL" > "$sfa"

  if [[ ! -s "$sfa" ]]; then
    echo "[SKIP] ${samp}: no spacers"
    continue
  fi

  # build DB once
  if [[ ! -f "${fa}.nhr" ]]; then
    echo "[mkdb] makeblastdb -in ${fa} -dbtype nucl"
    makeblastdb -in "$fa" -dbtype nucl >/dev/null
  fi

  echo "[MAP] ${samp}: spacers -> contigs"
  # best local hit per spacer on its own contigs
  blastn -db "$db" -query "$sfa" -task blastn-short -word_size 7 \
         -dust no -soft_masking false -evalue 1e-6 -max_target_seqs 5 \
         -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  | awk '
      {k=$1}
      !seen[k] || ($12+0) > best[k] { best[k]=$12; line[k]=$0; seen[k]=1 }  # keep highest bitscore
      END{ for (k in line) print line[k] }
    ' \
  | awk -v OFS="\t" '{print $1,$2}' > "$tsv"

  cat "$tsv" >> "$MAPALL"
done

# de-dup combined (prefer first)
awk -v OFS="\t" '!seen[$1]++' "$MAPALL" > "${MAPALL}.dedup"
mv -f "${MAPALL}.dedup" "$MAPALL"

echo "[OK] spacer→contig map -> $MAPALL"
head -n 5 "$MAPALL" || true
