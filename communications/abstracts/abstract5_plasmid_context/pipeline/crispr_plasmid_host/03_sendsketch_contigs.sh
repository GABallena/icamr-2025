#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

# Per-contig taxonomy with SendSketch (remote), with CHUNKING to avoid 100k limit.
# Inputs : filtered_contigs/*.fa
# Outputs: work/sendsketch/<sample>_contig_lineage.tsv (per-sample)
#          work/sendsketch/contig_lineage.all.tsv     (combined)
#
# Tunables (env):
#   CHUNK_NSEQ        number of contigs per chunk (default 40000; must be <100000)
#   SENDSKETCH_XMX    Java heap for sendsketch (default -Xmx4g)
#   RECORDS_PER_QUERY number of top hits to keep (default 1)

FILTERED_DIR="${FILTERED_DIR:-filtered_contigs}"
WORK="${WORK:-work}"
OUTDIR="${WORK}/sendsketch"
CHUNK_NSEQ="${CHUNK_NSEQ:-40000}"
SENDSKETCH_XMX="${SENDSKETCH_XMX:--Xmx4g}"
RECORDS_PER_QUERY="${RECORDS_PER_QUERY:-1}"

command -v sendsketch.sh >/dev/null || { echo "[ERR] sendsketch.sh not found in PATH"; exit 1; }
mkdir -p "${OUTDIR}"

shopt -s nullglob
samples=("${FILTERED_DIR}"/*.fa)
[ ${#samples[@]} -gt 0 ] || { echo "[ERR] No FASTA files in ${FILTERED_DIR}"; exit 1; }

COMBINED="${OUTDIR}/contig_lineage.all.tsv"
echo -e "sample\tcontig\tlineage" > "${COMBINED}"

parse_perseq () {
  # $1=infile(.txt or .gz)  $2=sample  $3=out_tsv (append mode if exists)
  local inpath="$1" sample="$2" outtsv="$3"
  python3 - "$inpath" "$sample" "$outtsv" << 'PY'
import sys, re, pathlib, gzip
src, sample, outp = sys.argv[1], sys.argv[2], sys.argv[3]
op = gzip.open if src.endswith('.gz') else open
q_re = re.compile(r'^Query:\s+(.+?)\s')

rows=[]; cur=None; header=False
with op(src, 'rt') as f:
    for line in f:
        if line.startswith('Query:'):
            m=q_re.match(line)
            cur = m.group(1) if m else None
            header=False
            continue
        if line.startswith('WKID'):  # header for a result block
            header=True
            continue
        if header and line.strip():
            parts=line.rstrip('\n').split('\t')
            tax = parts[-1] if len(parts)>1 else "NA"
            if cur: rows.append((cur, tax))
            header=False

mode = 'a' if pathlib.Path(outp).exists() else 'w'
with open(outp, mode) as o:
    if mode == 'w': o.write("contig\tlineage\n")
    for c,t in rows:
        o.write(f"{c}\t{t}\n")
print(f"[OK] parsed {len(rows)} contigs from {src} -> {outp}")
PY
}

for fa in "${samples[@]}"; do
  sample=$(basename "${fa}" .fa)
  chunkdir="${OUTDIR}/chunks/${sample}"
  mkdir -p "${chunkdir}"
  rm -f "${chunkdir}/"*.fa 2>/dev/null || true

  echo "[SPLIT] ${sample} → chunks of ${CHUNK_NSEQ} contigs"
  # Split FASTA into chunks with at most CHUNK_NSEQ sequences each
  awk -v max="${CHUNK_NSEQ}" -v prefix="${chunkdir}/${sample}.part" '
    BEGIN{RS=">"; ORS=""}
    NR>1{
      hdr=$0; sub(/\n.*/,"",hdr)
      seq=$0; sub(/^[^\n]*\n/,"",seq); gsub(/\n/,"",seq)
      n++; part=int((n-1)/max)+1
      file=sprintf("%s%03d.fa", prefix, part)
      print ">"hdr"\n"seq"\n" >> file
    }' "${fa}"

  contig_tsv="${OUTDIR}/${sample}_contig_lineage.tsv"
  rm -f "${contig_tsv}" 2>/dev/null || true

  shopt -s nullglob
  parts=( "${chunkdir}/${sample}.part"*.fa )
  if [ ${#parts[@]} -eq 0 ]; then
    echo "[WARN] No chunks produced for ${sample} (empty file?)"
    continue
  fi

  echo "[SKETCH] ${sample} (${#parts[@]} chunks)"
  for p in "${parts[@]}"; do
    basep=$(basename "${p}" .fa)
    outtxt="${OUTDIR}/${basep}_perseq.txt"

    # Per-contig remote SendSketch; keep top N records to minimize output
    # NOTE: Do NOT set refseq=t here (it’s the default server); do NOT exceed 100k sequences per request
    if ! sendsketch.sh ${SENDSKETCH_XMX} in="${p}" out="${outtxt}" mode=sequence persequence=t printall=t records="${RECORDS_PER_QUERY}" ow=t; then
      echo "[WARN] sendsketch failed for chunk ${basep}; skipping this chunk"
      continue
    fi

    # Some builds gzip output automatically; parse whichever exists
    if   [ -s "${outtxt}" ]; then parse_perseq "${outtxt}" "${sample}" "${contig_tsv}"
    elif [ -s "${outtxt}.gz" ]; then parse_perseq "${outtxt}.gz" "${sample}" "${contig_tsv}"
    else echo "[WARN] No per-contig output for chunk ${basep}"
    fi
  done

  if [ -s "${contig_tsv}" ]; then
    awk -v s="${sample}" 'BEGIN{FS=OFS="\t"} NR>1{print s,$1,$2}' "${contig_tsv}" >> "${COMBINED}"
  else
    echo "[WARN] No contig-level taxonomy for ${sample}"
  fi
done

echo "[OK] Combined contig map → ${COMBINED}"
head -n 5 "${COMBINED}" || true
