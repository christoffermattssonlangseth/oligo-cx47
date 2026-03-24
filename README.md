# Cx47 Oligo Notebook Scaffold

This project turns the question set in `Key questions, SM, V2.docx` into a notebook-first analysis scaffold for exploring `.h5ad` datasets related to oligodendroglial Cx47, demyelination, and remyelination.

The Word document is the scientific brief. The code here is the exploration layer:

- `config/dataset_catalog.csv`: candidate datasets pulled from the document
- `src/cx47_oligo/`: reusable helpers for `.h5ad` inspection, gene-panel matching, and question-specific scoring
- `src/cx47_oligo/scanpy_compat.py`: repo-local `scanpy` bootstrap with safe cache paths
- `src/cx47_oligo/scanpy_tools.py`: preprocessing, panel scoring, UMAP, and marker ranking helpers built on `scanpy`
- `src/cx47_oligo/exports.py`: helper functions for writing CSV, JSON, text, and figure outputs
- `notebooks/`: one inventory notebook plus one notebook per question
- `data/raw/`: place local `.h5ad` files here
- `scripts/import_catalog_datasets.py`: convert configured matrix-plus-metadata sources into `.h5ad`
- `results/`: tangible notebook outputs written to disk

## Questions Covered

1. Panglial network spatiotemporal remodeling during demyelination and repair
2. Cx47 and mitochondrial metabolic transition during remyelination
3. ER stress and reactive glial states during panglial remodeling
4. Transglial mitochondrial communication during neuroinflammation
5. Lipid metabolism and glymphatic or meningeal clearance pathways

## Setup

```bash
conda activate cellcharter
python -m pip install -e .
python scripts/import_catalog_datasets.py
jupyter lab
```

Then place one or more `.h5ad` files under `data/raw/` and open the notebooks in `notebooks/`.

For datasets listed in `config/dataset_catalog.csv` with `source_expression_path` and `source_metadata_path`, the importer can build `.h5ad` files directly from those tabular sources.

Each question notebook now writes concrete outputs under `results/<notebook>/<dataset>/`, including summary tables, grouped score tables, ranked markers, and figure files.

## Scanpy Runtime Note

Your `cellcharter` environment currently imports `scanpy` successfully only when `NUMBA_CACHE_DIR` and `MPLCONFIGDIR` point to writable directories. The project now handles this automatically through `cx47_oligo.scanpy_compat.import_scanpy()`, which places caches under `.cache/scanpy/` inside the repo before importing `scanpy`.

If you run standalone commands outside the notebooks, use the same pattern:

```python
from cx47_oligo.scanpy_compat import import_scanpy

sc = import_scanpy(".")
```

## Suggested Workflow

1. Start with `notebooks/00_dataset_inventory.ipynb`.
2. Inspect `.obs` metadata candidates for timepoint, condition, cell type, lesion state, and sample ID.
3. Run the question-specific notebook matching your current biological question.
4. Adjust the grouping columns in the notebook to the actual columns present in your dataset.

## Current Repo State

- The Word document is present.
- No `.h5ad` files are currently inside this project folder.
- The scaffold is ready to accept local `.h5ad` files when you add them.
