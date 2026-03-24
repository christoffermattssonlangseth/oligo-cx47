from __future__ import annotations

from array import array
from pathlib import Path
from typing import Any
import gzip

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def _open_text(path: str | Path):
    file_path = Path(path)
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt", encoding="utf-8", errors="replace")
    return open(file_path, "rt", encoding="utf-8", errors="replace")


def _make_unique(values: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    unique: list[str] = []
    for value in values:
        count = seen.get(value, 0)
        if count == 0:
            unique.append(value)
        else:
            unique.append(f"{value}-{count}")
        seen[value] = count + 1
    return unique


def read_metadata_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def import_expression_with_metadata(
    expression_path: str | Path,
    metadata_path: str | Path,
    dataset_id: str,
    title: str,
    value_type: str = "counts",
    zero_clip_abs: float = 0.0,
) -> ad.AnnData:
    expr_path = Path(expression_path)
    meta_path = Path(metadata_path)
    metadata = read_metadata_table(meta_path)

    with _open_text(expr_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        feature_label = header[0] or "gene"
        cell_ids = header[1:]
        if not cell_ids:
            raise ValueError(f"No cell identifiers found in header for {expr_path}")

        missing_meta = sorted(set(cell_ids) - set(metadata.index))
        if missing_meta:
            raise ValueError(
                f"Metadata is missing {len(missing_meta)} cell IDs from the expression matrix. "
                f"First missing ID: {missing_meta[0]}"
            )

        metadata = metadata.reindex(cell_ids).copy()
        metadata.index.name = "cell_id"
        metadata["source_cell_id"] = metadata.index

        indices = array("I")
        data = array("f")
        indptr = [0]
        genes: list[str] = []

        for line_number, line in enumerate(handle, start=2):
            gene, sep, values_text = line.partition("\t")
            if not sep:
                raise ValueError(f"Malformed expression row at line {line_number} in {expr_path}")

            values = np.fromstring(values_text, sep="\t", dtype=np.float32)
            if values.shape[0] != len(cell_ids):
                raise ValueError(
                    f"Line {line_number} in {expr_path} has {values.shape[0]} values but expected {len(cell_ids)}."
                )

            if zero_clip_abs > 0:
                values[np.abs(values) <= zero_clip_abs] = 0.0

            nz = np.flatnonzero(values)
            if nz.size:
                indices.frombytes(np.asarray(nz, dtype=np.uint32).tobytes())
                data.frombytes(np.asarray(values[nz], dtype=np.float32).tobytes())

            indptr.append(len(data))
            genes.append(gene)

    matrix = sparse.csc_matrix(
        (
            np.frombuffer(data, dtype=np.float32),
            np.frombuffer(indices, dtype=np.uint32),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=(len(cell_ids), len(genes)),
    )
    matrix.sort_indices()

    unique_var_names = _make_unique(genes)
    var = pd.DataFrame(
        {
            "gene_symbol": genes,
            "feature_type": feature_label,
        },
        index=pd.Index(unique_var_names, name="gene_id"),
    )

    adata = ad.AnnData(X=matrix, obs=metadata, var=var)
    adata.uns["source_dataset"] = {
        "dataset_id": dataset_id,
        "title": title,
        "expression_path": str(expr_path),
        "metadata_path": str(meta_path),
        "value_type": value_type,
        "zero_clip_abs": float(zero_clip_abs),
    }
    adata.uns["source_stats"] = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "nnz": int(matrix.nnz),
        "density": float(matrix.nnz / (adata.n_obs * adata.n_vars)),
    }
    return adata


def import_catalog_dataset(row: pd.Series, repo_root: str | Path) -> Path:
    required = ("dataset_id", "title", "local_path", "source_expression_path", "source_metadata_path")
    missing = [key for key in required if not row.get(key)]
    if missing:
        raise ValueError(f"Catalog row is missing required fields: {', '.join(missing)}")

    repo_path = Path(repo_root)
    output_path = repo_path / row["local_path"]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = import_expression_with_metadata(
        expression_path=row["source_expression_path"],
        metadata_path=row["source_metadata_path"],
        dataset_id=row["dataset_id"],
        title=row["title"],
        value_type=row.get("value_type", "counts") or "counts",
        zero_clip_abs=float(row.get("zero_clip_abs", 0.0) or 0.0),
    )
    adata.write_h5ad(output_path, compression="gzip")
    return output_path
