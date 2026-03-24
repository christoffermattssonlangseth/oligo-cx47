from __future__ import annotations

from pathlib import Path

from IPython.display import HTML
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .scanpy_compat import import_scanpy


PRIMARY = "#17465f"
SECONDARY = "#0f7a8a"
ACCENT = "#d88c32"
INK = "#11212d"


def configure_notebook_style(repo_root: str | Path | None = None):
    sc = import_scanpy(repo_root=repo_root)
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "#fbfcfd",
            "axes.edgecolor": "#d7dde2",
            "axes.labelcolor": INK,
            "axes.titleweight": "bold",
            "axes.titlesize": 15,
            "figure.titlesize": 16,
            "grid.color": "#dde4ea",
            "grid.linestyle": "-",
            "grid.alpha": 0.75,
            "xtick.color": INK,
            "ytick.color": INK,
            "font.size": 12,
        }
    )
    sc.set_figure_params(dpi=120, dpi_save=240, format="png", facecolor="white", frameon=False)
    return sc


def title_card_html(title: str, subtitle: str, dataset_hint: str = "") -> HTML:
    dataset_html = f"<p style='margin:10px 0 0 0;font-size:14px;opacity:0.9'>{dataset_hint}</p>" if dataset_hint else ""
    return HTML(
        f"""
        <div style="
            padding: 24px 28px;
            border-radius: 22px;
            background: linear-gradient(135deg, {PRIMARY} 0%, {SECONDARY} 58%, {ACCENT} 100%);
            color: white;
            box-shadow: 0 18px 40px rgba(23, 70, 95, 0.18);
            margin: 8px 0 18px 0;
        ">
            <div style="font-size: 26px; font-weight: 800; letter-spacing: 0.02em;">{title}</div>
            <p style="margin:10px 0 0 0;font-size:15px;line-height:1.5;max-width:880px;">{subtitle}</p>
            {dataset_html}
        </div>
        """
    )


def plot_group_counts(obs: pd.DataFrame, groupby: str, top_n: int = 25) -> pd.DataFrame:
    counts = obs[groupby].astype(str).value_counts().head(top_n).sort_values()
    fig, ax = plt.subplots(figsize=(8, max(4, 0.38 * len(counts))))
    counts.plot.barh(ax=ax, color=PRIMARY)
    ax.set_title(f"Cell counts by {groupby}")
    ax.set_xlabel("Cells")
    ax.set_ylabel("")
    sns.despine(ax=ax)
    plt.tight_layout()
    return counts.rename("n_cells").to_frame()


def plot_score_heatmap(obs: pd.DataFrame, groupby: str, score_columns: list[str]) -> pd.DataFrame:
    summary = obs.groupby(groupby, observed=False)[score_columns].mean().sort_index()
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(score_columns) + 2), max(4, 0.42 * len(summary) + 2)))
    sns.heatmap(summary, cmap="mako", center=0, linewidths=0.4, linecolor="white", ax=ax)
    ax.set_title(f"Mean panel scores by {groupby}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    return summary


def composition_table(
    obs: pd.DataFrame,
    row_group: str,
    col_group: str,
    top_n: int = 12,
    normalize: str = "row",
) -> pd.DataFrame:
    counts = pd.crosstab(obs[row_group].astype(str), obs[col_group].astype(str))
    if top_n and counts.shape[1] > top_n:
        keep = counts.sum(axis=0).sort_values(ascending=False).head(top_n).index
        other = counts.drop(columns=keep).sum(axis=1)
        counts = counts.loc[:, keep]
        if other.sum() > 0:
            counts["Other"] = other

    if normalize == "row":
        totals = counts.sum(axis=1).replace(0, np.nan)
        return counts.div(totals, axis=0).fillna(0)
    if normalize == "column":
        totals = counts.sum(axis=0).replace(0, np.nan)
        return counts.div(totals, axis=1).fillna(0)
    return counts


def plot_composition(obs: pd.DataFrame, row_group: str, col_group: str, top_n: int = 12) -> pd.DataFrame:
    fractions = composition_table(obs, row_group=row_group, col_group=col_group, top_n=top_n, normalize="row")
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(fractions.columns) + 4), max(4, 0.45 * len(fractions.index) + 2)))
    fractions.plot.bar(stacked=True, ax=ax, width=0.85, colormap="tab20")
    ax.set_title(f"{col_group} composition within {row_group}")
    ax.set_xlabel(row_group)
    ax.set_ylabel("Fraction")
    ax.legend(title=col_group, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    return fractions


def plot_score_distributions(obs: pd.DataFrame, groupby: str, score_columns: list[str], max_scores: int = 4) -> pd.DataFrame:
    selected_scores = list(score_columns[:max_scores])
    long_df = obs[[groupby, *selected_scores]].melt(id_vars=groupby, var_name="score", value_name="value")
    fig, ax = plt.subplots(figsize=(max(8, 2.4 * len(selected_scores)), 5))
    sns.boxplot(data=long_df, x="score", y="value", hue=groupby, ax=ax, fliersize=0.4)
    ax.set_title(f"Score distributions by {groupby}")
    ax.set_xlabel("")
    ax.set_ylabel("Score")
    ax.legend(title=groupby, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    return long_df


def plot_top_ranked_genes(rank_df: pd.DataFrame, top_n: int = 5) -> pd.DataFrame:
    score_col = "logfoldchanges" if "logfoldchanges" in rank_df.columns else "scores"
    top_df = rank_df.groupby("group", observed=False).head(top_n).copy()
    top_df["label"] = top_df["names"].astype(str)

    fig_height = max(4, 1.5 * top_df["group"].nunique())
    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.barplot(data=top_df, y="label", x=score_col, hue="group", dodge=False, ax=ax)
    ax.set_title(f"Top ranked genes by group ({score_col})")
    ax.set_xlabel(score_col)
    ax.set_ylabel("")
    ax.legend(title="group", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    return top_df


def plot_feature_heatmap(summary: pd.DataFrame, title: str, cmap: str = "mako", center: float | None = None) -> pd.DataFrame:
    plot_df = summary.copy()
    if isinstance(plot_df.index, pd.MultiIndex):
        plot_df.index = [" | ".join(map(str, idx)) for idx in plot_df.index]
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(plot_df.columns) + 2), max(4, 0.45 * len(plot_df.index) + 2)))
    sns.heatmap(plot_df, cmap=cmap, linewidths=0.4, linecolor="white", center=center, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    return summary


def plot_gene_contribution(summary: pd.DataFrame, title: str) -> pd.DataFrame:
    plot_df = summary.copy()
    if isinstance(plot_df.index, pd.MultiIndex):
        plot_df.index = [" | ".join(map(str, idx)) for idx in plot_df.index]

    fig, ax = plt.subplots(figsize=(max(8, 0.85 * len(plot_df.index) + 4), 5.4))
    plot_df.plot.bar(stacked=True, ax=ax, width=0.86, colormap="viridis")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of selected-gene signal")
    ax.set_ylim(0, 1.0)
    ax.legend(title="Gene", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    return summary


def plot_scatter_relationship(
    frame: pd.DataFrame,
    x: str,
    y: str,
    hue: str | None = None,
    style: str | None = None,
    title: str | None = None,
) -> pd.DataFrame:
    plot_df = frame.dropna(subset=[x, y]).copy()
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    sns.scatterplot(data=plot_df, x=x, y=y, hue=hue, style=style, ax=ax, s=42, alpha=0.8)
    ax.set_title(title or f"{y} vs {x}")
    sns.despine(ax=ax)
    plt.tight_layout()
    return plot_df
