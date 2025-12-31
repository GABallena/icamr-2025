#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail
WORK="${WORK:-work}"
MIN_SPACER_LEN="${MIN_SPACER_LEN:-15}"
mkdir -p "$WORK/spacers_fa"
SPALL="$WORK/spacers_fa/ALL.spacers.fasta"
> "$SPALL"

shopt -s nullglob
for rpt in "$WORK/minced/"*.report.txt; do
  base=$(basename "$rpt" .report.txt)
  sfa="$WORK/minced/${base}.spacers.fa"
  outnorm="$WORK/spacers_fa/${base}.normalized.fa"
  src="crt"
  if [ -s "$sfa" ] && grep -q '^>' "$sfa"; then src="sfa"; fi
  echo "[SPACERS] $base source=$src"

  python3 - "$src" "$rpt" "$sfa" "$outnorm" "$base" "$MIN_SPACER_LEN" << 'PY'
import sys,re
src, rpt, sfa, outfa, sample, minlen = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6])

def iter_sfa(path):
  header=None; seq=[]
  with open(path) as f:
    for line in f:
      line=line.strip()
      if not line: continue
      if line.startswith(">"):
        if header is not None:
          yield "".join(seq).upper()
        header=line[1:]; seq=[]
      else:
        seq.append(line)
    if header is not None:
      yield "".join(seq).upper()

def iter_crt(path):
  # match CRT table rows: " 175  REPEATSEQ  SPACERSEQ  [xx, yy]"
  row_re=re.compile(r'^\s*\d+\s+[ACGT\-]+\s+([ACGT]+)', re.I)
  with open(path) as f:
    for line in f:
      m=row_re.match(line)
      if m:
        yield m.group(1).upper()

def write_norm(seq_iter, out):
  i=0
  for s in seq_iter:
    s=re.sub(r'[^ACGTN]','',s)
    if len(s) < minlen: continue
    i+=1
    out.write(f">{sample}|NA|arr1|sp{i}|0-0\n{s}\n")

with open(outfa,"w") as o:
  it = iter_sfa(sfa) if src=="sfa" else iter_crt(rpt)
  write_norm(it, o)
PY

  cat "$outnorm" >> "$SPALL" || true
done

echo "[OK] Combined spacers → $SPALL"
grep -c '^>' "$SPALL" || true
