# Data Layout

Place local `.h5ad` files in `data/raw/`.

Suggested layout:

- `data/raw/`
- `data/processed/`
- `figures/`

The notebooks auto-scan `data/raw/` recursively for `.h5ad` files.

The dataset accession table in `config/dataset_catalog.csv` is only metadata. It does not download files for you.
