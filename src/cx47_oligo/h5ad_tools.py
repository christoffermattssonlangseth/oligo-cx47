from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .gene_panels import GENE_PANELS
from .questions import QUESTION_BANK

DEFAULT_DATASET_CATALOG = Path(__file__).resolve().parents[2] / "config" / "dataset_catalog.csv"
VAR_SYMBOL_COLUMNS = ("gene_symbol", "gene_symbols", "symbol", "features", "feature_name", "gene", "Gene")
OBS_KEYWORD_GROUPS = {
    "time_or_condition": ("time", "day", "week", "dpi", "condition", "stage", "lesion", "remyel", "demyel"),
    "cell_type_or_cluster": ("cell", "type", "cluster", "annotation", "annot", "label", "class", "subtype"),
    "sample_or_batch": ("sample", "batch", "replicate", "donor", "animal", "mouse", "patient", "library"),
    "spatial": ("region", "distance", "zone", "spatial", "x", "y", "slice", "section"),
}


def dataset_catalog(path: str | Path | None = None) -> pd.DataFrame:
    csv_path = Path(path) if path is not None else DEFAULT_DATASET_CATALOG
    return pd.read_csv(csv_path)


def discover_h5ad_files(root: str | Path) -> list[Path]:
    root_path = Path(root)
    if not root_path.exists():
        return []
    return sorted(root_path.rglob("*.h5ad"))


def load_h5ad(path: str | Path) -> ad.AnnData:
    file_path = Path(path)
    if not file_path.exists():
        raise FileNotFoundError(f"No file found at {file_path}")
    return ad.read_h5ad(file_path)


def adata_overview(adata: ad.AnnData) -> pd.DataFrame:
    rows = [
        ("n_obs", adata.n_obs),
        ("n_vars", adata.n_vars),
        ("obs_columns", len(adata.obs.columns)),
        ("var_columns", len(adata.var.columns)),
        ("layers", ", ".join(sorted(adata.layers.keys())) or "(none)"),
        ("obsm", ", ".join(sorted(adata.obsm.keys())) or "(none)"),
        ("uns_keys", ", ".join(sorted(adata.uns.keys())) or "(none)"),
        ("raw_present", adata.raw is not None),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def obs_column_summary(adata: ad.AnnData, max_examples: int = 4) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for column in adata.obs.columns:
        series = adata.obs[column]
        non_null = series.dropna()
        examples = ", ".join(non_null.astype(str).unique()[:max_examples])
        rows.append(
            {
                "column": column,
                "dtype": str(series.dtype),
                "n_unique": int(series.nunique(dropna=True)),
                "null_fraction": float(series.isna().mean()),
                "examples": examples,
            }
        )
    return pd.DataFrame(rows).sort_values(["n_unique", "column"], ascending=[True, True]).reset_index(drop=True)


def suggest_obs_columns(adata: ad.AnnData) -> pd.DataFrame:
    lowered = {column: column.lower() for column in adata.obs.columns}
    rows: list[dict[str, object]] = []
    for category, keywords in OBS_KEYWORD_GROUPS.items():
        matches = [column for column, lowered_name in lowered.items() if any(keyword in lowered_name for keyword in keywords)]
        rows.append({"category": category, "matches": ", ".join(matches) if matches else "(no match)"})
    return pd.DataFrame(rows)


def _gene_lookup(adata: ad.AnnData) -> dict[str, str]:
    mapping: dict[str, str] = {}
    var_names = [str(name) for name in adata.var_names]
    for var_name in var_names:
        mapping.setdefault(var_name.upper(), var_name)

    for column in VAR_SYMBOL_COLUMNS:
        if column not in adata.var.columns:
            continue
        aliases = adata.var[column].astype(str)
        for var_name, alias in zip(var_names, aliases):
            alias = alias.strip()
            if alias and alias.lower() != "nan":
                mapping.setdefault(alias.upper(), var_name)
    return mapping


def resolve_gene_symbols(adata: ad.AnnData, genes: Sequence[str]) -> tuple[list[str], list[str], dict[str, str]]:
    lookup = _gene_lookup(adata)
    present: list[str] = []
    missing: list[str] = []
    gene_to_var: dict[str, str] = {}
    for gene in genes:
        matched_var = lookup.get(gene.upper())
        if matched_var is None:
            missing.append(gene)
            continue
        present.append(gene)
        gene_to_var[gene] = matched_var
    return present, missing, gene_to_var


def panel_availability_table(adata: ad.AnnData, panel_names: Iterable[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for panel_name in panel_names:
        genes = GENE_PANELS[panel_name]
        present, missing, _ = resolve_gene_symbols(adata, genes)
        coverage = len(present) / len(genes) if genes else 0.0
        rows.append(
            {
                "panel": panel_name,
                "panel_size": len(genes),
                "present_count": len(present),
                "missing_count": len(missing),
                "coverage": round(coverage, 3),
                "present_genes": ", ".join(present),
                "missing_genes": ", ".join(missing),
            }
        )
    return pd.DataFrame(rows).sort_values(["coverage", "panel"], ascending=[False, True]).reset_index(drop=True)


def question_panel_table(adata: ad.AnnData, question_id: str) -> pd.DataFrame:
    question = QUESTION_BANK[question_id]
    return panel_availability_table(adata, question.panels)


def _matrix_to_numpy(matrix: object) -> np.ndarray:
    if sparse.issparse(matrix):
        values = matrix.toarray()
    else:
        values = np.asarray(matrix)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return values.astype(float, copy=False)


def _value_type(adata: ad.AnnData) -> str:
    source_dataset = adata.uns.get("source_dataset", {})
    return str(source_dataset.get("value_type", "")).lower()


def _cell_totals(adata: ad.AnnData, layer: str | None = None) -> np.ndarray:
    cache_key = "__cx47_total_counts" if layer is None else f"__cx47_total_counts_{layer}"
    if cache_key in adata.obs.columns:
        return adata.obs[cache_key].to_numpy(dtype=float, copy=False)

    matrix = adata.layers[layer] if layer is not None else adata.X
    if sparse.issparse(matrix):
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float, copy=False)
    else:
        totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float, copy=False)

    adata.obs[cache_key] = totals
    return totals


def expression_frame(
    adata: ad.AnnData,
    genes: Sequence[str],
    layer: str | None = None,
    normalize_counts: bool | str = "auto",
    target_sum: float = 1e4,
    log1p_normalized: bool = True,
) -> pd.DataFrame:
    present, _, gene_to_var = resolve_gene_symbols(adata, genes)
    if not present:
        return pd.DataFrame(index=adata.obs_names)

    ordered_vars = [gene_to_var[gene] for gene in present]
    view = adata[:, ordered_vars]
    matrix = view.layers[layer] if layer is not None else view.X
    values = _matrix_to_numpy(matrix)

    already_normalized = adata.uns.get("cx47_scanpy", {}).get("normalization") == "normalize_total+log1p"
    should_normalize = normalize_counts is True or (
        normalize_counts == "auto" and _value_type(adata) == "counts" and not already_normalized
    )
    if should_normalize:
        totals = _cell_totals(adata, layer=layer)
        scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        values = values * scale[:, None]
        if log1p_normalized:
            values = np.log1p(values)

    return pd.DataFrame(values, index=adata.obs_names, columns=present)


def _zscore_frame(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame.to_numpy(dtype=float, copy=True)
    means = np.nanmean(values, axis=0)
    stds = np.nanstd(values, axis=0)
    stds[stds == 0] = 1.0
    scaled = (values - means) / stds
    return pd.DataFrame(scaled, index=frame.index, columns=frame.columns)


def score_gene_panel(
    adata: ad.AnnData,
    panel_name: str,
    layer: str | None = None,
    standardize: bool = True,
    normalize_counts: bool | str = "auto",
) -> tuple[pd.Series, dict[str, object]]:
    genes = GENE_PANELS[panel_name]
    expr = expression_frame(adata, genes, layer=layer, normalize_counts=normalize_counts)
    if expr.empty:
        raise ValueError(f"No genes from panel '{panel_name}' were found in the AnnData object.")

    working = _zscore_frame(expr) if standardize and expr.shape[1] > 1 else expr
    score = working.mean(axis=1)
    score_name = f"{panel_name}_score"

    _, missing, _ = resolve_gene_symbols(adata, genes)
    metadata = {
        "panel": panel_name,
        "score_column": score_name,
        "genes_used": expr.shape[1],
        "genes_missing": len(missing),
        "missing_genes": ", ".join(missing),
        "normalization": "log1p_cpm" if normalize_counts is True or (normalize_counts == "auto" and _value_type(adata) == "counts") else "as_is",
    }
    return pd.Series(score, index=adata.obs_names, name=score_name), metadata


def score_question_panels(
    adata: ad.AnnData,
    question_id: str,
    layer: str | None = None,
    standardize: bool = True,
    normalize_counts: bool | str = "auto",
) -> dict[str, object]:
    question = QUESTION_BANK[question_id]
    obs_with_scores = adata.obs.copy()
    metadata_rows: list[dict[str, object]] = []
    score_columns: list[str] = []

    for panel_name in question.panels:
        try:
            score, metadata = score_gene_panel(
                adata,
                panel_name,
                layer=layer,
                standardize=standardize,
                normalize_counts=normalize_counts,
            )
            obs_with_scores[score.name] = score
            score_columns.append(score.name)
            metadata_rows.append(metadata)
        except ValueError:
            metadata_rows.append(
                {
                    "panel": panel_name,
                    "score_column": "(not scored)",
                    "genes_used": 0,
                    "genes_missing": len(GENE_PANELS[panel_name]),
                    "missing_genes": ", ".join(GENE_PANELS[panel_name]),
                }
            )

    return {
        "panel_summary": pd.DataFrame(metadata_rows),
        "obs_with_scores": obs_with_scores,
        "score_columns": score_columns,
    }


def mean_scores_by_group(obs_frame: pd.DataFrame, score_columns: Sequence[str], groupby: str) -> pd.DataFrame:
    if groupby not in obs_frame.columns:
        raise KeyError(f"'{groupby}' is not present in the provided observation frame.")
    if not score_columns:
        return pd.DataFrame()

    grouped = obs_frame.groupby(groupby, observed=False)[list(score_columns)].mean()
    counts = obs_frame.groupby(groupby, observed=False).size().rename("n_cells")
    return grouped.join(counts).sort_values("n_cells", ascending=False)
