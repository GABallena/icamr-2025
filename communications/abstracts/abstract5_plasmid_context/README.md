# Abstract 5 – Plasmid context / enrichment + CRISPR–plasmid host mapping (portfolio-safe)

Scripts:
- `abstract5_plasmid_enrichment.R` – enrichment-style summaries and basic modeling.

Pipeline helpers:
- `pipeline/crispr_plasmid_host/*.sh` – end-to-end or stepwise scripts to:
  1) extract spacers (PILER-CR or CRISPRCasFinder/minced-style reports),
  2) back-map spacers to contigs (optional),
  3) assign contig taxonomy (SendSketch),
  4) BLAST spacers vs a plasmid reference FASTA,
  5) join hits into a plasmid→host assignment table.

## Outputs
All scripts assume outputs go under a local `work/` directory (gitignored).
