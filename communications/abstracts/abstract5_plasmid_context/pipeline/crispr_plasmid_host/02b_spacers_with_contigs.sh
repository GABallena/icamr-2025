#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail
WORK="${WORK:-work}"
MIN_SPACER_LEN="${MIN_SPACER_LEN:-15}"
INRPT_DIR="${WORK}/minced"
OUTDIR="${WORK}/spacers_fa"
mkdir -p "${OUTDIR}"
SPALL="${OUTDIR}/ALL.spacers.fasta"
> "$SPALL"

shopt -s nullglob
for rpt in "${INRPT_DIR}"/*.report.txt; do
  base=$(basename "$rpt" .report.txt)
  outfa="${OUTDIR}/${base}.with_contigs.fa"
  echo "[SPACERS] $base"
  python3 - "$rpt" "$outfa" "$base" "$MIN_SPACER_LEN" << 'PY'
import sys, re
rpt, outfa, sample, minlen = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
seq_re = re.compile(r'(?i)^(?:Sequence(?:\s*[:])?|Seq(?:uence)?)\s+(\S+)')
row_re = re.compile(r'^\s*\d+\s+[ACGT\-]+\s+([ACGT]+)', re.I)
contig="NA"; spacers=[]
with open(rpt) as f:
    for line in f:
        m = seq_re.match(line)
        if m: contig = m.group(1); continue
        m = row_re.match(line)
        if m:
            s = re.sub(r'[^ACGTN]', '', m.group(1).upper())
            if len(s) >= minlen: spacers.append((contig, s))
with open(outfa,"w") as o:
    for i,(ctg,s) in enumerate(spacers,1):
        o.write(f">{sample}|{ctg}|arr1|sp{i}|0-0\n{s}\n")
print(f"[OK] {len(spacers)} -> {outfa}")
PY
  cat "$outfa" >> "$SPALL"
done

echo "[OK] Combined spacers → $SPALL"
grep -c '^>' "$SPALL" || true
