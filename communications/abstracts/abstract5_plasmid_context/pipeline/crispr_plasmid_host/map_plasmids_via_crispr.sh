#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail

RELABEL_DIR="relabeled_contigs"
PLASMIDS_FA="PLSDB_db/plasmids.slim.fasta"   # <-- point this to your plasmid DB
OUTDIR="work"
RAW="${OUTDIR}/crispr_raw"
SPC="${OUTDIR}/spacers_fa"
SSK="${OUTDIR}/sendsketch"
BLD="${OUTDIR}/blast"
RES="${OUTDIR}/results"

mkdir -p "$RAW" "$SPC" "$SSK" "$BLD" "$RES"

# 0) Build BLAST DB for plasmids (once)
if [ ! -f "${PLASMIDS_FA}.nin" ] && [ ! -f "${PLASMIDS_FA}.nhr" ]; then
  makeblastdb -in "$PLASMIDS_FA" -dbtype nucl
fi

# 1) For each contig fasta: detect CRISPR, extract spacers, assign taxonomy
for fa in ${RELABEL_DIR}/*_contigs.fa; do
  base=$(basename "$fa" .fa)

  # 1a) Detect CRISPR arrays (PILER-CR)
  # Output has arrays, repeats, and spacer sequences we will parse next
  pilercr -in "$fa" -out "${RAW}/${base}.pilercr.txt"

  # 1b) Extract spacers to FASTA with contig context
  python3 - "$fa" "${RAW}/${base}.pilercr.txt" > "${SPC}/${base}.spacers.fasta" << 'PY'
import sys,re
from pathlib import Path

# Args
fa_path = Path(sys.argv[1])          # contigs fasta (only for sample name)
pcr_path = Path(sys.argv[2])         # pilercr txt

sample = fa_path.stem
spacers = []
contig = None
array_id = 0
with open(pcr_path) as f:
    for line in f:
        line=line.strip()
        if line.startswith("Sequence "):
            # e.g., "Sequence SAMPLE-01_S1_contig_12345"
            contig = line.split("Sequence",1)[1].strip()
            array_id = 0
        elif line.startswith("Array"):
            array_id += 1
        elif line and re.match(r"^\d+\s+\d+\s+\d+\s+[ACGTNacgtn-]+", line):
            # pilercr spacer line: idx start end seq
            parts = line.split()
            idx, start, end, seq = parts[0], parts[1], parts[2], parts[3]
            # sanitize
            seq = seq.replace("-", "").upper()
            if len(seq) >= 15:  # keep short ones out; adjust if needed
                hdr = f">{sample}|{contig}|arr{array_id}|sp{idx}|{start}-{end}"
                spacers.append((hdr, seq))
# write
out = sys.stdout
for h,s in spacers:
    out.write(h + "\n")
    out.write(s + "\n")
PY

  # 1c) Assign taxonomy for contigs using BBMap SendSketch (k-mer taxon)
  # We only need taxonomy for contigs that carry arrays; SendSketch can run on full contigs
  # records=1 = best hit; add mode=tossjunk to reduce noise
  sendSketch_out="${SSK}/${base}.sendsketch.txt"
  sendsketch.sh in="$fa" out="$sendSketch_out" records=1

done

# 2) Concatenate all spacers into one file and BLAST to plasmids
cat ${SPC}/*.spacers.fasta > "${SPC}/ALL.spacers.fasta" || true

# BLAST: mirror paper’s permissive spacer matching (short reads)
# Paper used E=0.1 and filtered to <=2 mismatches and <=1 indel overall:{index=2}.
blastn -task blastn-short \
  -query "${SPC}/ALL.spacers.fasta" \
  -db "$PLASMIDS_FA" \
  -out "${BLD}/spacer_vs_plasmids.tsv" \
  -evalue 1e-1 \
  -word_size 7 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# 3) Parse SendSketch taxon → map back to spacer’s contig → join with BLAST hits
python3 - "${SSK}" "${BLD}/spacer_vs_plasmids.tsv" > "${RES}/plasmid_host_assignments.tsv" << 'PY'
import sys, re, glob
from collections import defaultdict

ssk_dir = sys.argv[1]
blast_tsv = sys.argv[2]

# Parse SendSketch (one best hit per contig) into dict: contig -> taxon_line
contig2tax = {}
for path in glob.glob(f"{ssk_dir}/*.sendsketch.txt"):
    with open(path) as f:
        for line in f:
            # Typical SendSketch table lines: name, taxID, name, ..., lineage...
            # We'll keep whole line as 'taxon' for transparency and extract genus later if needed.
            if line.startswith("#") or not line.strip(): continue
            parts = line.strip().split("\t")
            if len(parts) < 2: continue
            contig_name = parts[0]
            taxon = parts[-1]  # lineage string (last column)
            contig2tax[contig_name] = taxon

# Helper to get contig name from spacer header:
# >SAMPLE|CONTIG|arrX|spY|start-end
def parse_contig(qid):
    try:
        return qid.split("|")[1]
    except:
        return None

print("\t".join([
    "spacer_id","contig","host_lineage","plasmid_id","pident","aln_len",
    "mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"
]))

with open(blast_tsv) as f:
    for line in f:
        qid,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs = line.strip().split("\t")
        contig = parse_contig(qid)
        host = contig2tax.get(contig, "NA")
        print("\t".join([qid, contig or "NA", host, sid, pid, alen, mis, gap, qs, qe, ss, se, ev, bs]))
PY

echo "Done. See: ${RES}/plasmid_host_assignments.tsv"
