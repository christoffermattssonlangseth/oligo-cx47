from __future__ import annotations

from pathlib import Path
import re
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .gene_panels import GENE_PANELS
from .h5ad_tools import expression_frame, resolve_gene_symbols
from .questions import QUESTION_BANK
from .scanpy_compat import import_scanpy


def _value_type(adata: ad.AnnData) -> str:
    source_dataset = adata.uns.get("source_dataset", {})
    return str(source_dataset.get("value_type", "")).lower()


def ensure_scanpy_adata(
    adata: ad.AnnData,
    repo_root: str | Path | None = None,
    copy: bool = True,
    target_sum: float = 1e4,
    n_top_genes: int = 3000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    leiden_resolution: float = 0.5,
) -> ad.AnnData:
    sc = import_scanpy(repo_root=repo_root)
    work = adata.copy() if copy else adata
    work.var_names_make_unique()

    if "counts" not in work.layers:
        if sparse.issparse(work.X):
            work.layers["counts"] = work.X.copy()
        else:
            work.layers["counts"] = np.asarray(work.X).copy()

    value_type = _value_type(work)
    work.uns.setdefault("cx47_scanpy", {})
    work.uns["cx47_scanpy"]["input_value_type"] = value_type

    if value_type == "counts":
        sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
        work.uns["cx47_scanpy"]["normalization"] = "normalize_total+log1p"
    else:
        work.uns["cx47_scanpy"]["normalization"] = "as_is"

    if "highly_variable" not in work.var.columns:
        try:
            sc.pp.highly_variable_genes(work, n_top_genes=min(n_top_genes, work.n_vars), flavor="seurat")
        except Exception:
            work.var["highly_variable"] = True

    if "X_pca" not in work.obsm:
        sc.pp.pca(work, n_comps=min(n_pcs, max(2, work.n_vars - 1)), mask_var="highly_variable")
    if "neighbors" not in work.uns:
        sc.pp.neighbors(work, n_neighbors=n_neighbors, n_pcs=min(n_pcs, work.obsm["X_pca"].shape[1]))
    if "X_umap" not in work.obsm:
        sc.tl.umap(work)
    if "leiden" not in work.obs.columns:
        sc.tl.leiden(work, resolution=leiden_resolution, key_added="leiden")

    if work.raw is None:
        work.raw = work
    return work


def score_gene_panel_scanpy(
    adata: ad.AnnData,
    panel_name: str,
    repo_root: str | Path | None = None,
    score_name: str | None = None,
    use_raw: bool | None = None,
) -> dict[str, object]:
    sc = import_scanpy(repo_root=repo_root)
    genes = GENE_PANELS[panel_name]
    present, missing, gene_to_var = resolve_gene_symbols(adata, genes)
    if not present:
        raise ValueError(f"No genes from panel '{panel_name}' were found in the AnnData object.")

    score_column = score_name or f"{panel_name}_score"
    matched_var_names = [gene_to_var[gene] for gene in present]
    sc.tl.score_genes(adata, gene_list=matched_var_names, score_name=score_column, use_raw=use_raw, copy=False)
    return {
        "panel": panel_name,
        "score_column": score_column,
        "genes_used": len(present),
        "genes_missing": len(missing),
        "missing_genes": ", ".join(missing),
        "matched_var_names": ", ".join(matched_var_names),
        "method": "scanpy.tl.score_genes",
    }


def score_question_panels_scanpy(
    adata: ad.AnnData,
    question_id: str,
    repo_root: str | Path | None = None,
    use_raw: bool | None = None,
) -> dict[str, object]:
    question = QUESTION_BANK[question_id]
    metadata_rows: list[dict[str, object]] = []
    score_columns: list[str] = []

    for panel_name in question.panels:
        try:
            metadata = score_gene_panel_scanpy(
                adata,
                panel_name,
                repo_root=repo_root,
                score_name=f"{panel_name}_score",
                use_raw=use_raw,
            )
            metadata_rows.append(metadata)
            score_columns.append(metadata["score_column"])
        except ValueError:
            metadata_rows.append(
                {
                    "panel": panel_name,
                    "score_column": "(not scored)",
                    "genes_used": 0,
                    "genes_missing": len(GENE_PANELS[panel_name]),
                    "missing_genes": ", ".join(GENE_PANELS[panel_name]),
                    "method": "scanpy.tl.score_genes",
                }
            )

    return {"panel_summary": pd.DataFrame(metadata_rows), "score_columns": score_columns}


def question_marker_genes(question_id: str, top_per_panel: int = 5) -> list[str]:
    question = QUESTION_BANK[question_id]
    markers: list[str] = []
    seen: set[str] = set()
    for panel_name in question.panels:
        for gene in GENE_PANELS[panel_name][:top_per_panel]:
            if gene not in seen:
                seen.add(gene)
                markers.append(gene)
    return markers


def available_question_marker_genes(adata: ad.AnnData, question_id: str, top_per_panel: int = 5) -> list[str]:
    markers = question_marker_genes(question_id, top_per_panel=top_per_panel)
    present, _, gene_to_var = resolve_gene_symbols(adata, markers)
    return [gene_to_var[gene] for gene in present]


def rank_genes_groups_df(
    adata: ad.AnnData,
    groupby: str,
    repo_root: str | Path | None = None,
    method: str = "wilcoxon",
    key_added: str | None = None,
    n_genes: int = 25,
) -> pd.DataFrame:
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' is not present in adata.obs.")

    sc = import_scanpy(repo_root=repo_root)
    result_key = key_added or f"rank_genes_{groupby}"
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, key_added=result_key)
    return sc.get.rank_genes_groups_df(adata, group=None, key=result_key).head(n_genes * adata.obs[groupby].nunique())


def obs_projection_frame(
    adata: ad.AnnData,
    obs_columns: Sequence[str] | None = None,
    score_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    selected_columns: list[str] = []
    for column in list(obs_columns or []) + list(score_columns or []):
        if column in adata.obs.columns and column not in selected_columns:
            selected_columns.append(column)

    frame = adata.obs[selected_columns].copy() if selected_columns else pd.DataFrame(index=adata.obs_names)
    frame.index.name = "cell_id"

    if "X_umap" in adata.obsm:
        frame["umap_1"] = adata.obsm["X_umap"][:, 0]
        frame["umap_2"] = adata.obsm["X_umap"][:, 1]

    if "leiden" in adata.obs.columns and "leiden" not in frame.columns:
        frame["leiden"] = adata.obs["leiden"].astype(str)

    return frame


def obs_keyword_mask(adata: ad.AnnData, column: str, keywords: Sequence[str]) -> pd.Series:
    if column not in adata.obs.columns:
        raise KeyError(f"'{column}' is not present in adata.obs.")
    pattern = "|".join(re.escape(keyword) for keyword in keywords)
    return adata.obs[column].astype(str).str.contains(pattern, case=False, regex=True, na=False)


def assign_keyword_families(values: pd.Series, keyword_map: dict[str, Sequence[str]], default: str = "Other") -> pd.Series:
    families = pd.Series(default, index=values.index, dtype="object")
    values_as_text = values.astype(str)
    for family, keywords in keyword_map.items():
        pattern = "|".join(re.escape(keyword) for keyword in keywords)
        mask = values_as_text.str.contains(pattern, case=False, regex=True, na=False)
        families.loc[mask] = family
    return families


def grouped_gene_expression(
    adata: ad.AnnData,
    groupby: str,
    genes: Sequence[str],
    normalize_counts: bool | str = "auto",
) -> pd.DataFrame:
    return grouped_gene_expression_stats(
        adata,
        groupby=groupby,
        genes=genes,
        normalize_counts=normalize_counts,
    )["mean_expression"]


def grouped_gene_expression_stats(
    adata: ad.AnnData,
    groupby: str,
    genes: Sequence[str],
    normalize_counts: bool | str = "auto",
    expr_threshold: float = 0.0,
) -> dict[str, pd.DataFrame]:
    if groupby not in adata.obs.columns:
        raise KeyError(f"'{groupby}' is not present in adata.obs.")

    expr = expression_frame(adata, genes, normalize_counts=normalize_counts)
    if expr.empty:
        empty = pd.DataFrame()
        return {
            "n_cells": empty,
            "mean_expression": empty,
            "pct_expressing": empty,
            "gene_contribution": empty,
        }

    group_labels = adata.obs[groupby].astype(str)
    gene_columns = list(expr.columns)

    expr_with_groups = expr.copy()
    expr_with_groups[groupby] = group_labels.values

    mean_expression = expr_with_groups.groupby(groupby, observed=False)[gene_columns].mean()

    detected = expr.gt(expr_threshold).astype(float)
    detected[groupby] = group_labels.values
    pct_expressing = detected.groupby(groupby, observed=False)[gene_columns].mean().mul(100.0)

    counts = expr_with_groups.groupby(groupby, observed=False).size().rename("n_cells").to_frame()
    order = counts["n_cells"].sort_values(ascending=False).index

    mean_expression = mean_expression.reindex(order)
    pct_expressing = pct_expressing.reindex(order)
    counts = counts.reindex(order)

    row_totals = mean_expression.sum(axis=1).replace(0, np.nan)
    gene_contribution = mean_expression.div(row_totals, axis=0).fillna(0.0)

    return {
        "n_cells": counts,
        "mean_expression": mean_expression,
        "pct_expressing": pct_expressing,
        "gene_contribution": gene_contribution,
    }
