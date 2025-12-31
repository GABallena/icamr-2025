# CRISPR–plasmid host mapping helpers (scripts only)

This folder contains helper scripts for a CRISPR-spacer-based plasmid→host linkage workflow.

There are two main approaches:
1) **PILER-CR** (`pilercr_hosts.sh`) – run PILER-CR on contigs, extract spacers, BLAST spacers vs plasmids, join with taxonomy.
2) **CRISPRCasFinder** (`ccf_plasmid_hosts.sh`) – similar, but uses CRISPRCasFinder outputs.

For more control / resumability, run the stepwise scripts:
- `02_build_spacers.sh` (or `02b_spacers_with_contigs.sh`) → build normalized spacer FASTA(s)
- `02b_backmap_spacers_to_contigs.sh` → optional: map spacers back to sample contigs
- `03_sendsketch_contigs.sh` → contig taxonomy
- `04_blast_spacers_resume.sh` → resumable BLAST against plasmid FASTA
- `05_join_results_contig.sh` or `05_join_results_contig_with_map.sh` → final joined table

## Expected working directories
By default, scripts write to:
- `${WORK:-work}/...`

Keep `work/` out of git.

## Dependencies
Common:
- bash, python3
- BLAST+ (`makeblastdb`, `blastn`)
- BBMap `sendsketch.sh`

Depending on approach:
- `pilercr` (PILER-CR)
- `CRISPRCasFinder.pl` + dependencies (perl, RNAfold, hmmsearch, prodigal)
