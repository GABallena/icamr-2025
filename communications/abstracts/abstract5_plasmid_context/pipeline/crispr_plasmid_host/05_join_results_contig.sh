#!/usr/bin/env bash
# Portfolio-safe copy: paths/IDs generalized; databases/outputs not included.
#
set -euo pipefail
WORK="${WORK:-work}"
BLTSV="${WORK}/blast/spacer_vs_plasmids.tsv"
MAPCTG="${WORK}/sendsketch/contig_lineage.all.tsv"
OUT="${WORK}/results/plasmid_host_assignments.contig.tsv"
MM_MAX="${MM_MAX:-2}"
GAPOPEN_MAX="${GAPOPEN_MAX:-1}"
mkdir -p "${WORK}/results"

python3 - "$BLTSV" "$MAPCTG" "$OUT" "$MM_MAX" "$GAPOPEN_MAX" << 'PY'
import sys, csv
blt, mapf, outf, mmx, gox = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])

# sample+contig → lineage
sc2line = {}
with open(mapf) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        sc2line[(row['sample'], row['contig'])] = row['lineage']

def parse_spacer_id(qid):
    parts = qid.split('|')
    sample = parts[0] if len(parts)>0 else "NA"
    contig = parts[1] if len(parts)>1 else "NA"
    return sample, contig

with open(outf,"w") as o, open(blt) as f:
    o.write("\t".join(["spacer_id","sample","contig","host_lineage",
                       "plasmid_id","pident","aln_len","mismatch","gapopen",
                       "qstart","qend","sstart","send","evalue","bitscore"])+"\n")
    for line in f:
        cols=line.rstrip("\n").split("\t")
        if len(cols)<12: continue
        qid,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs = cols
        if int(float(mis))>mmx or int(float(gap))>gox: 
            continue
        s,c = parse_spacer_id(qid)
        host = sc2line.get((s,c), "NA")
        o.write("\t".join([qid,s,c,host,sid,pid,alen,mis,gap,qs,qe,ss,se,ev,bs])+"\n")
print(f"[OK] wrote {outf}")
PY

echo "[OK] Contig-level results → ${OUT}"
echo "[TIP] unique plasmid–host pairs:"
cut -f4,5 "${OUT}" | sort -u | wc -l
