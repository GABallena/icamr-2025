# icamr-2025

Materials, figures, and analysis scripts for an AMR-focused conference submission.

This repository is intended to be **portfolio-safe**:
- no raw sequencing reads
- no sensitive identifiers (sample IDs, exact sites, institutions)
- figures and summaries should be generated from **public or anonymized** inputs only

## What’s inside (typical)
- `scripts/` — data wrangling, stats, and figure generation
- `figures/` — rendered figures for posters/slides (exported outputs)
- `notes/` — short writeups, figure legends, and reproducibility notes

## How to run
1. Create a local (untracked) `data/` folder and place your input tables there.
2. Run scripts from the repo root so relative paths work.

Example:
```bash
Rscript scripts/01_build_tables.R
Rscript scripts/02_make_figures.R
```

## License
MIT — see `LICENSE`.
