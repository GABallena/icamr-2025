#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail
WORK="${WORK:-work}"
BLTSV="${BLTSV:-${WORK}/blast/spacer_vs_plasmids.tsv}"
MAPCTG="${WORK}/sendsketch/contig_lineage.all.tsv"
SP2CTG="${WORK}/spacers_fa/spacer_to_contig.tsv"
OUT="${WORK}/results/plasmid_host_assignments.contig.tsv"
MM_MAX="${MM_MAX:-2}"
GAPOPEN_MAX="${GAPOPEN_MAX:-1}"
mkdir -p "${WORK}/results"

python3 - <<'PY' "$BLTSV" "$MAPCTG" "$SP2CTG" "$OUT" "$MM_MAX" "$GAPOPEN_MAX"
import sys, csv, re
blt, mapf, spmap, outf, mmx, gox = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6])

# spacer -> contig (from back-map)
sp2c = {}
with open(spmap) as f:
    for row in csv.reader(f, delimiter="\t"):
        if row: sp2c[row[0]] = row[1]

# (sample, contig) -> lineage
sc2 = {}
with open(mapf) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        sc2[(row['sample'], row['contig'])] = row['lineage']

def parse_qid(qid):
    parts = qid.split('|')
    sample = parts[0] if len(parts) > 0 else "NA"
    contig = parts[1] if len(parts) > 1 else "NA"
    return sample, contig

bad=0; kept=0
with open(outf, "w") as out, open(blt) as f:
    out.write("\t".join(["spacer_id","sample","contig","host_lineage",
                         "plasmid_id","pident","aln_len","mismatch","gapopen",
                         "qstart","qend","sstart","send","evalue","bitscore"])+"\n")
    for line in f:
        line=line.rstrip("\n")
        if not line: continue
        # robust split: tabs first; if not 12 cols, split on any whitespace
        cols = re.split(r"\t+", line)
        if len(cols) != 12:
            cols = line.split()
        if len(cols) != 12:
            bad += 1
            continue

        qid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = cols

        if int(float(mismatch)) > mmx or int(float(gapopen)) > gox:
            continue

        sample, contig = parse_qid(qid)
        if contig == "NA":
            contig = sp2c.get(qid, "NA")

        host = sc2.get((sample, contig), "NA")

        out.write("\t".join([qid, sample, contig, host, sseqid, pident, length, mismatch, gapopen,
                             qstart, qend, sstart, send, evalue, bitscore]) + "\n")
        kept += 1

print(f"[OK] wrote {outf} (kept={kept}, skipped_bad_rows={bad})")
PY

echo "[OK] Contig-level results → ${OUT}"
echo "[TIP] unique plasmid–host pairs:"
cut -f4,5 "${OUT}" | sort -u | wc -l
