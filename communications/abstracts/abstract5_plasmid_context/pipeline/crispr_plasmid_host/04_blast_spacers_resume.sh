#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

# Resumable BLAST of spacers vs PLSDB nucleotide FASTA (chunked)
WORK="${WORK:-work}"
Q="${WORK}/spacers_fa/ALL.spacers.fasta"
DB="${PLASMIDS_DB:-PLSDB_db/plasmids.slim.fasta}"
BLAST_DIR="${WORK}/blast"
OUT_ALL="${BLAST_DIR}/spacer_vs_plasmids.tsv"

# Tunables
CHUNK_NSEQ="${CHUNK_NSEQ:-50000}"     # spacers per chunk
NUM_THREADS="${NUM_THREADS:-4}"       # threads per chunk
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-10}"
EVALUE="${EVALUE:-1e-5}"

mkdir -p "${BLAST_DIR}"

[ -s "$Q" ] || { echo "[ERR] No spacers: $Q is empty"; exit 1; }
command -v blastn >/dev/null || { echo "[ERR] blastn not found"; exit 1; }

# Ensure nucleotide DB is formatted
if [[ ! -f "${DB}.nhr" ]]; then
  echo "[mkdb] makeblastdb -in ${DB} -dbtype nucl"
  makeblastdb -in "${DB}" -dbtype nucl
fi

# Split queries into chunks (resumable: keep if already present)
CHUNK_DIR="${BLAST_DIR}/chunks"
mkdir -p "${CHUNK_DIR}"
if ! compgen -G "${CHUNK_DIR}/spacers.part"*.fa >/dev/null; then
  echo "[SPLIT] $(basename "$Q") → chunks of ${CHUNK_NSEQ} seqs"
  awk -v max="${CHUNK_NSEQ}" -v prefix="${CHUNK_DIR}/spacers.part" '
    BEGIN{RS=">"; ORS=""}
    NR>1{
      hdr=$0; sub(/\n.*/,"",hdr)
      seq=$0; sub(/^[^\n]*\n/,"",seq)
      n++; part=int((n-1)/max)+1
      file=sprintf("%s%03d.fa", prefix, part)
      print ">"hdr"\n"seq >> file
    }' "$Q"
fi

mapfile -t parts < <(ls "${CHUNK_DIR}/spacers.part"*.fa | sort)
echo "[INFO] ${#parts[@]} chunks to process"

for p in "${parts[@]}"; do
  base=$(basename "$p" .fa)
  out="${BLAST_DIR}/${base}.tsv"
  tmp="${out}.tmp"

  if [[ -s "$out" ]]; then
    echo "[SKIP] ${base} (exists)"
    continue
  fi

  echo "[BLAST] ${base}"
  stdbuf -oL blastn \
    -db "${DB}" \
    -query "${p}" \
    -task blastn-short \
    -word_size 7 \
    -dust no \
    -soft_masking false \
    -evalue "${EVALUE}" \
    -max_target_seqs "${MAX_TARGET_SEQS}" \
    -num_threads "${NUM_THREADS}" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    > "${tmp}"

  mv -f "${tmp}" "${out}"
done

# Merge all chunk TSVs
echo "[MERGE] → ${OUT_ALL}"
: > "${OUT_ALL}"
for t in "${BLAST_DIR}"/spacers.part*.tsv; do
  [ -s "$t" ] && cat "$t" >> "${OUT_ALL}"
done

echo "[OK] ${OUT_ALL}"
wc -l "${OUT_ALL}" || true
head -n 3 "${OUT_ALL}" || true
