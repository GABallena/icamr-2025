#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

# -------- CONFIG --------
RELABEL_DIR="${RELABEL_DIR:-relabeled_contigs}"       # input contigs
OUTBASE="${OUTBASE:-CCF_OUT}"                          # CRISPRCasFinder outputs
WORK="${WORK:-work}"                                   # scratch + results
PLASMIDS_FA="${PLASMIDS_FA:-PLSDB_db/plasmids.slim.fasta}"  # PLSDB FASTA (existing)

# Stringency similar to paper’s spacer matching
EVALUE="${EVALUE:-1e-1}"
WORDSIZE="${WORDSIZE:-7}"
MM_MAX="${MM_MAX:-2}"        # <= mismatches
GAPOPEN_MAX="${GAPOPEN_MAX:-1}"  # <= gap opens
MIN_SPACER_LEN="${MIN_SPACER_LEN:-15}"

# -------- PREP --------
mkdir -p "${OUTBASE}" "${WORK}"/{spacers_fa,sendsketch,blast,results}

# Check binaries
need() { command -v "$1" >/dev/null 2>&1 || { echo "[ERR] $1 not found"; exit 1; }; }
need perl
need blastn
need sendsketch.sh
need RNAfold
need hmmsearch
need prodigal

# CRISPRCasFinder script path (assumed in PATH; adjust if needed)
CCF_BIN="${CCF_BIN:-CRISPRCasFinder.pl}"
if ! command -v "${CCF_BIN}" >/dev/null 2>&1; then
  echo "[ERR] CRISPRCasFinder.pl not found."
  echo "      Install it and ensure it's on PATH, or set CCF_BIN=/path/to/CRISPRCasFinder.pl"
  exit 1
fi
fi

# BLAST DB for PLSDB
if [ ! -f "${PLASMIDS_FA}.nhr" ] && [ ! -f "${PLASMIDS_FA}.00.nhr" ]; then
  echo "[INFO] Building BLAST DB for ${PLASMIDS_FA}"
  makeblastdb -in "${PLASMIDS_FA}" -dbtype nucl
fi

# -------- 1) RUN CRISPRCasFinder PER SAMPLE --------
echo "[INFO] Running CRISPRCasFinder on contigs in ${RELABEL_DIR}"
shopt -s nullglob
for fa in "${RELABEL_DIR}"/*_contigs.fa; do
  base=$(basename "$fa" .fa)
  outdir="${OUTBASE}/${base}"
  if [ -d "${outdir}" ] && [ -f "${outdir}/CRISPR-Cas_summary.tsv" ]; then
    echo "[SKIP] ${base} appears done"
  else
    echo "[RUN] ${base}"
    perl "${CCF_BIN}" -in "$fa" -meta -keep -cas -def General -out "${outdir}"
  fi
done

# -------- 2) GATHER SPACERS AND NORMALIZE HEADERS --------
# We’ll create a spacer FASTA with headers: >SAMPLE|CONTIG|arrN|spM|start-end
SPALL="${WORK}/spacers_fa/ALL.spacers.fasta"
> "${SPALL}"

# CRISPRCasFinder usually emits spacer FASTAs named like *_Spacers*.fa in result dirs.
echo "[INFO] Collecting spacers"
for sfa in $(find "${OUTBASE}" -type f -iname '*Spacers*.fa'); do
  sample=$(basename "$(dirname "$sfa")")
  # Normalize headers using a flexible Python parser that tries to pull contig + indices from header lines
  python3 - << 'PY' "$sfa" "$sample" "${MIN_SPACER_LEN}"
import sys, re
from pathlib import Path
from itertools import groupby

sfa   = Path(sys.argv[1])
sample= sys.argv[2]
minlen= int(sys.argv[3])

def read_fasta(p):
    with open(p) as f:
        faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
        for header in faiter:
            h = next(header).strip()[1:]
            seq = "".join(s.strip() for s in next(faiter))
            yield h, seq

# Heuristics to extract contig/array/spacer/start-end from CCF headers.
# We accept anything like:
#   spacer_#|SequenceID|start-end ...
#   spacer_#_arrN_... seqid=..., start=..., end=...
# We’ll try multiple patterns; fallback to sample as contig if missing.
for h,seq in read_fasta(sfa):
    if len(seq) < minlen: 
        continue
    contig = None; arr="arr1"; sp="sp1"; pos="0-0"
    # try patterns
    m = re.search(r'(?:\||\s)([^|\s]+)\|(\d+)-(\d+)', h)  # ... |CONTIG|start-end
    if m:
        contig = m.group(1); pos=f"{m.group(2)}-{m.group(3)}"
    else:
        m = re.search(r'seqid=([^,\s]+).*?start=(\d+).*?end=(\d+)', h)
        if m:
            contig = m.group(1); pos=f"{m.group(2)}-{m.group(3)}"
        else:
            # last resort: try to grab a token that looks like a contig id
            m = re.search(r'(scaffold|contig)[^|\s,]*', h)
            if m: contig = m.group(0)
    ma = re.search(r'(?:arr|array)[_:=]?(\d+)', h, re.I)
    if ma: arr=f"arr{ma.group(1)}"
    ms = re.search(r'(?:spacer|sp)[_:=]?(\d+)', h, re.I)
    if ms: sp=f"sp{ms.group(1)}"
    if not contig:
        contig = sample  # fallback
    print(f">{sample}|{contig}|{arr}|{sp}|{pos}")
    # ensure sequence uppercase and no gaps
    print(re.sub(r'[^ACGTN]', '', seq.upper()))
PY
done > "${SPALL}"

# Safety
if [ ! -s "${SPALL}" ]; then
  echo "[WARN] No spacers collected. Check CRISPRCasFinder outputs."
fi

# -------- 3) TAXONOMY (SendSketch) FOR CONTIG FASTAS --------
echo "[INFO] Running SendSketch taxonomy"
for fa in "${RELABEL_DIR}"/*_contigs.fa; do
  base=$(basename "$fa" .fa)
  out="\"${WORK}/sendsketch/${base}.txt\""
  # If already exists with content, skip
  if [ -s ${WORK}/sendsketch/${base}.txt ]; then
    echo "[SKIP] taxonomy ${base}"
  else
    sendsketch.sh in="$fa" out="${WORK}/sendsketch/${base}.txt" records=1
  fi
done

# -------- 4) BLAST spacers → PLSDB --------
echo "[INFO] BLAST spacers to ${PLASMIDS_FA}"
BLTSV="${WORK}/blast/spacer_vs_plasmids.tsv"
blastn -task blastn-short \
  -query "${SPALL}" \
  -db "${PLASMIDS_FA}" \
  -out "${BLTSV}" \
  -evalue "${EVALUE}" -word_size "${WORDSIZE}" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# -------- 5) JOIN: qseqid -> contig -> taxonomy, then filter by mismatches/gaps --------
echo "[INFO] Joining hits with taxonomy"
python3 - << 'PY' "${WORK}/sendsketch" "${BLTSV}" "${WORK}/results/plasmid_host_assignments.tsv" "${MM_MAX}" "${GAPOPEN_MAX}"
import sys, glob

ssk_dir = sys.argv[1]
blast_tsv = sys.argv[2]
out_tsv = sys.argv[3]
mm_max = int(sys.argv[4])
gap_max = int(sys.argv[5])

# load sendsketch (contig -> lineage)
contig2tax = {}
for path in glob.glob(f"{ssk_dir}/*.txt"):
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.strip().split("\t")
            if len(parts) < 2: continue
            contig = parts[0]
            lineage = parts[-1]
            contig2tax[contig] = lineage

def parse_qid(qid):
    # >SAMPLE|CONTIG|arrN|spM|start-end
    parts = qid.split("|")
    if len(parts) >= 2:
        return parts[0], parts[1]
    return "NA","NA"

with open(out_tsv, "w") as out, open(blast_tsv) as f:
    out.write("\t".join(["spacer_id","sample","contig","host_lineage","plasmid_id",
                         "pident","aln_len","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])+"\n")
    for line in f:
        qid,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs = line.strip().split("\t")
        mis_i = int(float(mis)); gap_i = int(float(gap))
        if mis_i > mm_max or gap_i > gap_max:
            continue
        sample, contig = parse_qid(qid)
        host = contig2tax.get(contig, "NA")
        out.write("\t".join([qid,sample,contig,host,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs])+"\n")
PY

echo "[OK] Results: ${WORK}/results/plasmid_host_assignments.tsv"
echo "[TIP] Count unique plasmid→host links:"
echo "cut -f4,5 ${WORK}/results/plasmid_host_assignments.tsv | sort -u | wc -l"
