#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

# ---------- CONFIG (override via env) ----------
RELABEL_DIR="${RELABEL_DIR:-relabeled_contigs}"              # input contigs (*.fa)
WORK="${WORK:-work}"                                         # scratch + outputs
PLASMIDS_FA="${PLASMIDS_FA:-PLSDB_db/plasmids.slim.fasta}"   # your PLSDB FASTA
EVALUE="${EVALUE:-1e-1}"
WORDSIZE="${WORDSIZE:-7}"
MM_MAX="${MM_MAX:-2}"            # <= mismatches
GAPOPEN_MAX="${GAPOPEN_MAX:-1}"  # <= gap opens
MIN_SPACER_LEN="${MIN_SPACER_LEN:-15}"

mkdir -p "${WORK}"/{crispr_raw,spacers_fa,sendsketch,blast,results}

need(){ command -v "$1" >/dev/null 2>&1 || { echo "[ERR] $1 not found"; exit 1; }; }
need pilercr
need sendsketch.sh
need blastn
need python3

# Build BLAST DB if missing
if [ ! -f "${PLASMIDS_FA}.nhr" ] && [ ! -f "${PLASMIDS_FA}.00.nhr" ]; then
  echo "[INFO] Building BLAST DB for ${PLASMIDS_FA}"
  makeblastdb -in "${PLASMIDS_FA}" -dbtype nucl
fi

SPALL="${WORK}/spacers_fa/ALL.spacers.fasta"
> "${SPALL}"

shopt -s nullglob
echo "[INFO] Running PILER-CR on contigs in ${RELABEL_DIR}"
for fa in "${RELABEL_DIR}"/*_contigs.fa; do
  base=$(basename "$fa" .fa)
  outtxt="${WORK}/crispr_raw/${base}.pilercr.txt"
  if [ ! -s "${outtxt}" ]; then
    # PILER-CR writes to -out text file (not FASTA)
    pilercr -in "$fa" -out "${outtxt}" -noinfo -quiet -minarray 3 -mincons 0.9 -minrepeat 16 -maxrepeat 64 -minspacer 8 -maxspacer 64 -minrepeatratio 0.9 -minspacerratio 0.75
  fi

  # Extract spacers from PILER-CR text → FASTA with normalized headers
  python3 - "$fa" "${outtxt}" "${MIN_SPACER_LEN}" >> "${SPALL}" << 'PY'
import sys,re
from pathlib import Path

fa, ptxt, minlen = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
sample = fa.stem  # e.g., SAMPLE-01_S1_contigs
contig = None
arr_id = 0

with open(ptxt) as f:
    for ln in f:
        ln = ln.strip()
        if ln.startswith("Sequence "):
            contig = ln.split("Sequence",1)[1].strip()
            arr_id = 0
        elif ln.startswith("Array"):
            arr_id += 1
        elif ln and re.match(r"^\d+\s+\d+\s+\d+\s+[ACGTNacgtn-]+", ln):
            idx, start, end, seq = ln.split()[:4]
            seq = seq.replace("-", "").upper()
            if len(seq) >= minlen:
                hdr = f">{sample}|{contig}|arr{arr_id}|sp{idx}|{start}-{end}"
                print(hdr); print(seq)
PY
done

if [ ! -s "${SPALL}" ]; then
  echo "[WARN] No spacers found. Exiting."; exit 0
fi

# Taxonomy (SendSketch) per contig file
echo "[INFO] SendSketch taxonomy"
for fa in "${RELABEL_DIR}"/*_contigs.fa; do
  base=$(basename "$fa" .fa)
  out="\"${WORK}/sendsketch/${base}.txt\""
  [ -s ${WORK}/sendsketch/${base}.txt ] || sendsketch.sh in="$fa" out="${WORK}/sendsketch/${base}.txt" records=1
done

# BLAST spacers → plasmids
echo "[INFO] BLAST spacers vs ${PLASMIDS_FA}"
BLTSV="${WORK}/blast/spacer_vs_plasmids.tsv"
blastn -task blastn-short \
  -query "${SPALL}" -db "${PLASMIDS_FA}" -out "${BLTSV}" \
  -evalue "${EVALUE}" -word_size "${WORDSIZE}" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# Join with taxonomy + filter by mismatches & gaps
echo "[INFO] Join hits with host taxonomy"
python3 - << 'PY' "${WORK}/sendsketch" "${BLTSV}" "${WORK}/results/plasmid_host_assignments.tsv" "${MM_MAX}" "${GAPOPEN_MAX}"
import sys,glob
ssk,blt,out,mmx,gox = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

# contig -> lineage from SendSketch tables
contig2tax={}
for p in glob.glob(f"{ssk}/*.txt"):
    with open(p) as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip(): continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            contig, lineage = parts[0], parts[-1]
            contig2tax[contig] = lineage

def parse_qid(q):
    # >SAMPLE|CONTIG|arrN|spM|start-end
    pr = q.split("|")
    return (pr[0], pr[1]) if len(pr)>=2 else ("NA","NA")

with open(out,"w") as o, open(blt) as f:
    o.write("\t".join(["spacer_id","sample","contig","host_lineage","plasmid_id",
                       "pident","aln_len","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])+"\n")
    for line in f:
        qid,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs = line.strip().split("\t")
        if int(float(mis))>mmx or int(float(gap))>gox: continue
        sample, contig = parse_qid(qid)
        host = contig2tax.get(contig, "NA")
        o.write("\t".join([qid,sample,contig,host,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs])+"\n")
PY

echo "[OK] Results: ${WORK}/results/plasmid_host_assignments.tsv"
echo "[TIP] Unique plasmid→host combos:"
echo "cut -f4,5 ${WORK}/results/plasmid_host_assignments.tsv | sort -u | wc -l"
