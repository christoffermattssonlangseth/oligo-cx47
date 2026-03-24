from __future__ import annotations

import inspect
import sys
from pathlib import Path

import nbformat as nbf


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
NOTEBOOKS_DIR = REPO_ROOT / "notebooks"
sys.path.insert(0, str(SRC_ROOT))

from cx47_oligo.questions import ORDERED_QUESTION_IDS, QUESTION_BANK


DEFAULT_NOTEBOOK_SETTINGS = {
    "q1": {
        "dataset_path": 'REPO_ROOT / "data" / "raw" / "jakel_et_al.h5ad"',
        "group_columns": '["Lesion", "Celltypes"]',
        "dataset_label": "jakel_et_al",
    },
    "q2": {
        "dataset_path": 'REPO_ROOT / "data" / "raw" / "jakel_et_al.h5ad"',
        "group_columns": '["Lesion", "Celltypes"]',
        "dataset_label": "jakel_et_al",
    },
    "q3": {
        "dataset_path": 'REPO_ROOT / "data" / "raw" / "falcao_et_al.h5ad"',
        "group_columns": '["Group", "Renamed_clusternames"]',
        "dataset_label": "falcao_et_al",
    },
    "q4": {
        "dataset_path": 'REPO_ROOT / "data" / "raw" / "falcao_et_al.h5ad"',
        "group_columns": '["Group", "Renamed_clusternames"]',
        "dataset_label": "falcao_et_al",
    },
    "q5": {
        "dataset_path": 'REPO_ROOT / "data" / "raw" / "jakel_et_al.h5ad"',
        "group_columns": '["Lesion", "Celltypes"]',
        "dataset_label": "jakel_et_al",
    },
}


IMPORTS = """from pathlib import Path
import sys

from IPython.display import display

CANDIDATES = [Path.cwd().resolve(), Path.cwd().resolve().parent]
REPO_ROOT = next((path for path in CANDIDATES if (path / "src").exists() and (path / "config").exists()), None)
if REPO_ROOT is None:
    raise RuntimeError("Could not locate the repository root. Start Jupyter from the repo root or notebooks directory.")

if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

for module_name in [name for name in list(sys.modules) if name == "cx47_oligo" or name.startswith("cx47_oligo.")]:
    del sys.modules[module_name]

from cx47_oligo.gene_panels import GENE_PANELS
from cx47_oligo.exports import ensure_results_dir, save_current_figure, save_json, save_table, save_text, slugify
from cx47_oligo.h5ad_tools import (
    adata_overview,
    dataset_catalog,
    discover_h5ad_files,
    expression_frame,
    load_h5ad,
    mean_scores_by_group,
    obs_column_summary,
    question_panel_table,
    suggest_obs_columns,
)
from cx47_oligo.notebook_viz import (
    composition_table,
    configure_notebook_style,
    plot_composition,
    plot_feature_heatmap,
    plot_gene_contribution,
    plot_group_counts,
    plot_scatter_relationship,
    plot_score_distributions,
    plot_score_heatmap,
    plot_top_ranked_genes,
    title_card_html,
)
from cx47_oligo.questions import QUESTION_BANK
from cx47_oligo.scanpy_tools import (
    available_question_marker_genes,
    assign_keyword_families,
    ensure_scanpy_adata,
    grouped_gene_expression,
    grouped_gene_expression_stats,
    obs_projection_frame,
    obs_keyword_mask,
    rank_genes_groups_df,
    score_question_panels_scanpy,
)

sc = configure_notebook_style(REPO_ROOT)
"""


def markdown_cell(text: str) -> nbf.NotebookNode:
    return nbf.v4.new_markdown_cell(inspect.cleandoc(text) + "\n")


def code_cell(code: str) -> nbf.NotebookNode:
    return nbf.v4.new_code_cell(inspect.cleandoc(code) + "\n")


def build_inventory_notebook() -> nbf.NotebookNode:
    nb = nbf.v4.new_notebook()
    nb.cells = [
        markdown_cell(
            """
            # Dataset Inventory

            Use this notebook to audit the local `.h5ad` files you have placed in `data/raw/`, inspect the dataset catalog extracted from the Word document, and choose the dataset you want to analyze next.
            """
        ),
        code_cell(IMPORTS),
        code_cell(
            """
            display(
                title_card_html(
                    "Dataset Inventory",
                    "Audit the imported h5ad objects, inspect their metadata structure, and choose the strongest starting point for each biological question.",
                    dataset_hint="Current default datasets: Jäkel for MS lesion-state questions, Falcão for EAE stress and reactive-glia questions.",
                )
            )
            """
        ),
        code_cell(
            """
            catalog = dataset_catalog()
            display(catalog)
            """
        ),
        code_cell(
            """
            DATA_ROOT = REPO_ROOT / "data" / "raw"
            LOCAL_H5AD_FILES = discover_h5ad_files(DATA_ROOT)
            LOCAL_H5AD_FILES
            """
        ),
        code_cell(
            """
            DATASET_PATH = LOCAL_H5AD_FILES[0] if LOCAL_H5AD_FILES else None
            DATASET_PATH
            """
        ),
        code_cell(
            """
            RESULTS_DIR = ensure_results_dir(REPO_ROOT, "00_dataset_inventory", DATASET_PATH)
            FIGURES_DIR = RESULTS_DIR / "figures"
            FIGURES_DIR.mkdir(parents=True, exist_ok=True)
            sc.settings.figdir = FIGURES_DIR
            print("Results will be written to:", RESULTS_DIR)
            save_table(catalog, RESULTS_DIR / "dataset_catalog_snapshot.csv", index=False)
            save_json({"dataset_path": str(DATASET_PATH) if DATASET_PATH else "", "results_dir": str(RESULTS_DIR)}, RESULTS_DIR / "run_metadata.json")
            """
        ),
        code_cell(
            """
            adata = load_h5ad(DATASET_PATH) if DATASET_PATH else None
            if adata is None:
                print("No .h5ad file found under data/raw yet.")
            else:
                overview = adata_overview(adata)
                suggestions = suggest_obs_columns(adata)
                obs_summary = obs_column_summary(adata)
                display(overview)
                display(suggestions)
                display(obs_summary.head(25))
                save_table(overview, RESULTS_DIR / "adata_overview.csv", index=False)
                save_table(suggestions, RESULTS_DIR / "suggested_obs_columns.csv", index=False)
                save_table(obs_summary, RESULTS_DIR / "obs_column_summary.csv", index=False)
                print("scanpy version:", sc.__version__)
            """
        ),
        code_cell(
            """
            DEFAULT_GROUP = next(
                (column for column in ["Celltypes", "Renamed_clusternames", "Lesion", "Group", "Condition"] if adata is not None and column in adata.obs.columns),
                None,
            )

            if adata is not None and DEFAULT_GROUP is not None:
                group_counts = plot_group_counts(adata.obs, DEFAULT_GROUP)
                display(group_counts)
                save_table(group_counts, RESULTS_DIR / f"cell_counts_by_{slugify(DEFAULT_GROUP)}.csv")
                save_current_figure(FIGURES_DIR / f"cell_counts_by_{slugify(DEFAULT_GROUP)}.png")
            """
        ),
        markdown_cell(
            """
            ## Next Step

            Once you identify a useful dataset and the relevant `.obs` columns for timepoint, condition, and cell state, move to one of the question-specific notebooks in this folder.
            """
        ),
    ]
    return nb


def build_question_specific_cells(question_id: str) -> list[nbf.NotebookNode]:
    if question_id == "q1":
        return [
            markdown_cell(
                """
                ## Q1 Targeted Views

                These cells focus on oligodendroglial state composition and connexin behavior across lesion states.
                """
            ),
            code_cell(
                """
                OLIGO_COLUMN = next((column for column in ["Celltypes", "Renamed_clusternames", "leiden"] if adata_sc is not None and column in adata_sc.obs.columns), None)
                adata_oligo = None

                if adata_sc is not None and OLIGO_COLUMN is not None:
                    oligo_mask = obs_keyword_mask(adata_sc, OLIGO_COLUMN, ["oligo", "opc", "cop", "imolg", "nfol", "mol"])
                    adata_oligo = adata_sc[oligo_mask].copy()
                    print("Oligodendroglial cells:", adata_oligo.n_obs)

                if adata_oligo is not None and adata_oligo.n_obs > 0 and PRIMARY_GROUP is not None:
                    conn_summary = grouped_gene_expression(adata_oligo, PRIMARY_GROUP, ["GJC2", "GJB1", "GJA1"])
                    if not conn_summary.empty:
                        display(conn_summary)
                        save_table(conn_summary, RESULTS_DIR / f"q1_oligo_connexin_means_by_{slugify(PRIMARY_GROUP)}.csv")
                        display(plot_feature_heatmap(conn_summary, title=f"Oligodendroglial connexins by {PRIMARY_GROUP}", cmap="crest"))
                        save_current_figure(FIGURES_DIR / f"q1_oligo_connexin_heatmap_by_{slugify(PRIMARY_GROUP)}.png")

                        present_connexins = [gene for gene in ["GJC2", "GJB1", "GJA1"] if gene in conn_summary.columns]
                        if present_connexins:
                            sc.pl.dotplot(
                                adata_oligo,
                                var_names=present_connexins,
                                groupby=PRIMARY_GROUP,
                                standard_scale="var",
                                dendrogram=False,
                                save=f"_q1_oligo_connexins_{slugify(PRIMARY_GROUP)}.png",
                            )

                    if SECONDARY_GROUP is not None and SECONDARY_GROUP in adata_oligo.obs.columns:
                        oligo_comp = plot_composition(adata_oligo.obs, PRIMARY_GROUP, SECONDARY_GROUP)
                        display(oligo_comp)
                        save_table(oligo_comp, RESULTS_DIR / f"q1_oligo_composition_{slugify(SECONDARY_GROUP)}_within_{slugify(PRIMARY_GROUP)}.csv")
                        save_current_figure(FIGURES_DIR / f"q1_oligo_composition_{slugify(SECONDARY_GROUP)}_within_{slugify(PRIMARY_GROUP)}.png")
                else:
                    print("Need oligodendroglial annotations and a valid grouping column for Q1-specific plots.")
                """
            ),
        ]

    if question_id == "q2":
        return [
            markdown_cell(
                """
                ## Q2 Targeted Views

                These cells focus on the relationship between Cx47-linked signals and mitochondrial programs within oligodendroglia.
                """
            ),
            code_cell(
                """
                OLIGO_COLUMN = next((column for column in ["Celltypes", "Renamed_clusternames", "leiden"] if adata_sc is not None and column in adata_sc.obs.columns), None)
                adata_oligo = None

                if adata_sc is not None and OLIGO_COLUMN is not None:
                    oligo_mask = obs_keyword_mask(adata_sc, OLIGO_COLUMN, ["oligo", "opc", "cop", "imolg", "nfol", "mol"])
                    adata_oligo = adata_sc[oligo_mask].copy()
                    print("Oligodendroglial cells:", adata_oligo.n_obs)

                if adata_oligo is not None and adata_oligo.n_obs > 0 and PRIMARY_GROUP is not None:
                    q2_frame = adata_oligo.obs[[column for column in [PRIMARY_GROUP, "panglial_connexins_score", "mitochondrial_oxphos_score", "mitochondrial_biogenesis_score", "metabolic_coupling_score"] if column in adata_oligo.obs.columns]].copy()
                    gjc2_expr = expression_frame(adata_oligo, ["GJC2"])
                    if not gjc2_expr.empty:
                        q2_frame["GJC2_expr"] = gjc2_expr.iloc[:, 0].values

                    summary_cols = [column for column in ["GJC2_expr", "panglial_connexins_score", "mitochondrial_oxphos_score", "mitochondrial_biogenesis_score", "metabolic_coupling_score"] if column in q2_frame.columns]
                    q2_summary = q2_frame.groupby(PRIMARY_GROUP, observed=False)[summary_cols].mean()
                    display(q2_summary)
                    save_table(q2_summary, RESULTS_DIR / f"q2_cx47_mito_means_by_{slugify(PRIMARY_GROUP)}.csv")
                    display(plot_feature_heatmap(q2_summary, title=f"Cx47 and mitochondrial programs by {PRIMARY_GROUP}", cmap="rocket"))
                    save_current_figure(FIGURES_DIR / f"q2_cx47_mito_heatmap_by_{slugify(PRIMARY_GROUP)}.png")

                    scatter_y = next((column for column in ["mitochondrial_oxphos_score", "mitochondrial_biogenesis_score"] if column in q2_frame.columns), None)
                    scatter_x = "GJC2_expr" if "GJC2_expr" in q2_frame.columns else ("panglial_connexins_score" if "panglial_connexins_score" in q2_frame.columns else None)
                    if scatter_x is not None and scatter_y is not None:
                        scatter_df = plot_scatter_relationship(
                            q2_frame,
                            x=scatter_x,
                            y=scatter_y,
                            hue=PRIMARY_GROUP,
                            title=f"{scatter_y} vs {scatter_x} in oligodendroglia",
                        )
                        display(scatter_df.head())
                        save_table(scatter_df, RESULTS_DIR / f"q2_scatter_{slugify(scatter_x)}_vs_{slugify(scatter_y)}.csv", index=False)
                        save_current_figure(FIGURES_DIR / f"q2_scatter_{slugify(scatter_x)}_vs_{slugify(scatter_y)}.png")
                else:
                    print("Need oligodendroglial annotations and a valid grouping column for Q2-specific plots.")
                """
            ),
        ]

    if question_id == "q3":
        return [
            markdown_cell(
                """
                ## Q3 Targeted Views

                These cells focus on ER stress, inflammatory programs, and reactive glial-family structure in the Falcão-style EAE setting.
                """
            ),
            code_cell(
                """
                Q3_CELL_COLUMN = next((column for column in ["Renamed_clusternames", "Celltypes", "leiden"] if adata_sc is not None and column in adata_sc.obs.columns), None)

                if adata_sc is not None and Q3_CELL_COLUMN is not None and "Group" in adata_sc.obs.columns:
                    q3_frame = adata_sc.obs[[Q3_CELL_COLUMN, "Group", *[column for column in ["er_stress_upr_score", "inflammatory_cytokines_score", "microglial_activation_score", "panglial_connexins_score"] if column in adata_sc.obs.columns]]].copy()
                    q3_frame["cell_family"] = assign_keyword_families(
                        q3_frame[Q3_CELL_COLUMN],
                        {
                            "Microglia": ["migl", "micro"],
                            "Oligodendroglia": ["oligo", "opc", "cop", "imolg", "nfol", "mol"],
                            "Astrocytes": ["astro"],
                            "Vascular/VLMC": ["vlmc", "vascular", "endo", "pericyte"],
                        },
                    )

                    relevant = [column for column in ["er_stress_upr_score", "inflammatory_cytokines_score", "microglial_activation_score", "panglial_connexins_score"] if column in q3_frame.columns]
                    q3_summary = q3_frame.groupby(["Group", "cell_family"], observed=False)[relevant].mean()
                    display(q3_summary)
                    save_table(q3_summary, RESULTS_DIR / "q3_group_family_score_summary.csv")
                    display(plot_feature_heatmap(q3_summary, title="ER stress and reactive glia by Group and family", cmap="mako"))
                    save_current_figure(FIGURES_DIR / "q3_group_family_score_heatmap.png")

                    family_comp = plot_composition(q3_frame, "Group", "cell_family")
                    display(family_comp)
                    save_table(family_comp, RESULTS_DIR / "q3_cell_family_composition_by_group.csv")
                    save_current_figure(FIGURES_DIR / "q3_cell_family_composition_by_group.png")

                    if "er_stress_upr_score" in q3_frame.columns and "microglial_activation_score" in q3_frame.columns:
                        scatter_df = plot_scatter_relationship(
                            q3_frame,
                            x="er_stress_upr_score",
                            y="microglial_activation_score",
                            hue="Group",
                            style="cell_family",
                            title="ER stress vs microglial activation",
                        )
                        display(scatter_df.head())
                        save_table(scatter_df, RESULTS_DIR / "q3_er_stress_vs_microglial_activation.csv", index=False)
                        save_current_figure(FIGURES_DIR / "q3_er_stress_vs_microglial_activation.png")
                else:
                    print("Need a cluster-like annotation column and Group metadata for Q3-specific plots.")
                """
            ),
        ]

    return []


def build_question_notebook(question_id: str) -> nbf.NotebookNode:
    question = QUESTION_BANK[question_id]
    panel_list = ", ".join(question.panels)
    obs_keywords = ", ".join(question.recommended_obs_keywords)
    dataset_ids = ", ".join(question.dataset_ids)
    optional_label = "Optional question." if question.optional else "Core question."
    defaults = DEFAULT_NOTEBOOK_SETTINGS[question_id]

    nb = nbf.v4.new_notebook()
    nb.cells = [
        markdown_cell(
            f"""
            # {question.identifier.upper()}. {question.short_title}

            **Question**: {question.full_question}

            **Why this notebook exists**: {question.rationale}

            **Relevant panels**: `{panel_list}`

            **Recommended `.obs` keywords to look for**: `{obs_keywords}`

            **Suggested datasets from the brief**: `{dataset_ids}`

            **Scope**: {optional_label}
            """
        ),
        code_cell(IMPORTS),
        code_cell(
            f"""
            display(
                title_card_html(
                    "{question.identifier.upper()}. {question.short_title}",
                    "{question.full_question}",
                    dataset_hint="Suggested default dataset: {defaults['dataset_label']}",
                )
            )
            """
        ),
        markdown_cell(
            """
            ## Dataset And Metadata
            """
        ),
        code_cell(
            """
            DATA_ROOT = REPO_ROOT / "data" / "raw"
            LOCAL_H5AD_FILES = discover_h5ad_files(DATA_ROOT)
            LOCAL_H5AD_FILES
            """
        ),
        code_cell(
            f"""
            DATASET_PATH = {defaults["dataset_path"]}
            if not DATASET_PATH.exists():
                DATASET_PATH = LOCAL_H5AD_FILES[0] if LOCAL_H5AD_FILES else None
            DATASET_PATH
            """
        ),
        code_cell(
            f"""
            RESULTS_DIR = ensure_results_dir(REPO_ROOT, "{question.notebook_slug}", DATASET_PATH)
            FIGURES_DIR = RESULTS_DIR / "figures"
            FIGURES_DIR.mkdir(parents=True, exist_ok=True)
            sc.settings.figdir = FIGURES_DIR
            print("Results will be written to:", RESULTS_DIR)
            save_json(
                {{
                    "question_id": "{question.identifier}",
                    "question_title": "{question.short_title}",
                    "dataset_path": str(DATASET_PATH) if DATASET_PATH else "",
                    "group_columns": {defaults["group_columns"]},
                    "results_dir": str(RESULTS_DIR),
                }},
                RESULTS_DIR / "run_metadata.json",
            )
            """
        ),
        code_cell(
            """
            adata = load_h5ad(DATASET_PATH) if DATASET_PATH else None
            if adata is None:
                print("No .h5ad file found under data/raw yet. Add one or point DATASET_PATH to a local file.")
            else:
                overview = adata_overview(adata)
                display(overview)
                save_table(overview, RESULTS_DIR / "adata_overview.csv", index=False)
                print("scanpy version:", sc.__version__)
            """
        ),
        code_cell(
            f"""
            if adata is not None:
                suggested = suggest_obs_columns(adata)
                obs_summary = obs_column_summary(adata)
                panel_table = question_panel_table(adata, "{question_id}")
                display(suggested)
                display(obs_summary.head(30))
                display(panel_table)
                save_table(suggested, RESULTS_DIR / "suggested_obs_columns.csv", index=False)
                save_table(obs_summary, RESULTS_DIR / "obs_column_summary.csv", index=False)
                save_table(panel_table, RESULTS_DIR / "question_panel_availability.csv", index=False)
            """
        ),
        markdown_cell(
            """
            ## Scanpy Preparation
            """
        ),
        code_cell(
            f"""
            GROUP_COLUMNS = {defaults["group_columns"]}
            adata_sc = None
            results = None

            if adata is not None:
                adata_sc = ensure_scanpy_adata(adata, repo_root=REPO_ROOT, copy=True)
                results = score_question_panels_scanpy(adata_sc, "{question_id}", repo_root=REPO_ROOT)
                display(results["panel_summary"])
                if results["score_columns"]:
                    display(adata_sc.obs[results["score_columns"]].head())
                save_table(results["panel_summary"], RESULTS_DIR / "panel_score_summary.csv", index=False)
                projection = obs_projection_frame(adata_sc, obs_columns=GROUP_COLUMNS, score_columns=results["score_columns"])
                save_table(projection, RESULTS_DIR / "cell_level_scores_umap.csv")
                save_text(
                    "\\n".join(
                        [
                            f"Question: {question.full_question}",
                            f"Dataset: {{DATASET_PATH}}",
                            f"Results directory: {{RESULTS_DIR}}",
                            f"Group columns: {{', '.join(GROUP_COLUMNS)}}",
                            f"Score columns: {{', '.join(results['score_columns'])}}",
                        ]
                    ),
                    RESULTS_DIR / "analysis_summary.txt",
                )
            """
        ),
        markdown_cell(
            """
            ## Score Structure
            """
        ),
        code_cell(
            """
            if adata_sc is not None and results is not None:
                for groupby in GROUP_COLUMNS:
                    if groupby in adata_sc.obs.columns:
                        print(f"Score means by {groupby}")
                        grouped_scores = mean_scores_by_group(adata_sc.obs, results["score_columns"], groupby)
                        display(grouped_scores)
                        save_table(grouped_scores, RESULTS_DIR / f"score_means_by_{slugify(groupby)}.csv")
            """
        ),
        code_cell(
            """
            PRIMARY_GROUP = next((column for column in GROUP_COLUMNS if adata_sc is not None and column in adata_sc.obs.columns), None)
            SECONDARY_GROUP = next((column for column in GROUP_COLUMNS[1:] if adata_sc is not None and column in adata_sc.obs.columns), None)

            if adata_sc is not None and results is not None and PRIMARY_GROUP is not None:
                group_counts = plot_group_counts(adata_sc.obs, PRIMARY_GROUP)
                display(group_counts)
                save_table(group_counts, RESULTS_DIR / f"cell_counts_by_{slugify(PRIMARY_GROUP)}.csv")
                save_current_figure(FIGURES_DIR / f"cell_counts_by_{slugify(PRIMARY_GROUP)}.png")

                score_heatmap = plot_score_heatmap(adata_sc.obs, PRIMARY_GROUP, results["score_columns"])
                display(score_heatmap)
                save_table(score_heatmap, RESULTS_DIR / f"score_heatmap_values_by_{slugify(PRIMARY_GROUP)}.csv")
                save_current_figure(FIGURES_DIR / f"score_heatmap_by_{slugify(PRIMARY_GROUP)}.png")

                score_dist = plot_score_distributions(adata_sc.obs, PRIMARY_GROUP, results["score_columns"])
                display(score_dist.head())
                save_table(score_dist, RESULTS_DIR / f"score_distributions_by_{slugify(PRIMARY_GROUP)}.csv", index=False)
                save_current_figure(FIGURES_DIR / f"score_distributions_by_{slugify(PRIMARY_GROUP)}.png")

                if SECONDARY_GROUP is not None:
                    composition = plot_composition(adata_sc.obs, PRIMARY_GROUP, SECONDARY_GROUP)
                    display(composition)
                    save_table(composition, RESULTS_DIR / f"composition_{slugify(SECONDARY_GROUP)}_within_{slugify(PRIMARY_GROUP)}.csv")
                    save_current_figure(FIGURES_DIR / f"composition_{slugify(SECONDARY_GROUP)}_within_{slugify(PRIMARY_GROUP)}.png")
                    composition_counts = composition_table(adata_sc.obs, PRIMARY_GROUP, SECONDARY_GROUP, normalize="none")
                    save_table(composition_counts, RESULTS_DIR / f"composition_counts_{slugify(SECONDARY_GROUP)}_within_{slugify(PRIMARY_GROUP)}.csv")
            else:
                print("Set GROUP_COLUMNS so at least one valid grouping column is available.")
            """
        ),
        markdown_cell(
            """
            ## Marker-Level Views
            """
        ),
        code_cell(
            f"""
            PANEL_TO_INSPECT = QUESTION_BANK["{question_id}"].panels[0]
            MARKER_GENES = available_question_marker_genes(adata_sc, "{question_id}", top_per_panel=3) if adata_sc is not None else []

            if adata is not None:
                panel_expr = expression_frame(adata, GENE_PANELS[PANEL_TO_INSPECT])
                print(f"Matched {{panel_expr.shape[1]}} genes in {{PANEL_TO_INSPECT}}")
                display(panel_expr.head())
                save_table(panel_expr.head(100), RESULTS_DIR / f"panel_expression_preview_{{PANEL_TO_INSPECT}}.csv")
                save_text("\\n".join(MARKER_GENES), RESULTS_DIR / "marker_genes_used.txt")
            """
        ),
        markdown_cell(
            """
            ## Cell-Type-Specific Gene Expression
            """
        ),
        code_cell(
            """
            GENE_BREAKDOWN_GROUP = next(
                (
                    column
                    for column in ["Renamed_clusternames", PRIMARY_GROUP, SECONDARY_GROUP]
                    if column is not None and adata_sc is not None and column in adata_sc.obs.columns
                ),
                None,
            )
            GENE_BREAKDOWN_GENES = MARKER_GENES[: min(8, len(MARKER_GENES))]

            if adata_sc is not None and GENE_BREAKDOWN_GROUP is not None and GENE_BREAKDOWN_GENES:
                gene_breakdown = grouped_gene_expression_stats(
                    adata_sc,
                    groupby=GENE_BREAKDOWN_GROUP,
                    genes=GENE_BREAKDOWN_GENES,
                )

                display(gene_breakdown["n_cells"])
                save_table(gene_breakdown["n_cells"], RESULTS_DIR / f"gene_breakdown_cell_counts_by_{slugify(GENE_BREAKDOWN_GROUP)}.csv")
                save_text("\\n".join(GENE_BREAKDOWN_GENES), RESULTS_DIR / f"gene_breakdown_genes_by_{slugify(GENE_BREAKDOWN_GROUP)}.txt")

                mean_expression = gene_breakdown["mean_expression"]
                pct_expressing = gene_breakdown["pct_expressing"]
                gene_contribution = gene_breakdown["gene_contribution"]

                display(mean_expression)
                save_table(mean_expression, RESULTS_DIR / f"marker_mean_expression_by_{slugify(GENE_BREAKDOWN_GROUP)}.csv")
                display(plot_feature_heatmap(mean_expression, title=f"Mean marker expression by {GENE_BREAKDOWN_GROUP}", cmap="crest"))
                save_current_figure(FIGURES_DIR / f"marker_mean_expression_by_{slugify(GENE_BREAKDOWN_GROUP)}.png")

                display(pct_expressing)
                save_table(pct_expressing, RESULTS_DIR / f"marker_pct_expressing_by_{slugify(GENE_BREAKDOWN_GROUP)}.csv")
                display(plot_feature_heatmap(pct_expressing, title=f"Percent of cells expressing markers by {GENE_BREAKDOWN_GROUP}", cmap="flare"))
                save_current_figure(FIGURES_DIR / f"marker_pct_expressing_by_{slugify(GENE_BREAKDOWN_GROUP)}.png")

                display(gene_contribution)
                save_table(gene_contribution, RESULTS_DIR / f"marker_gene_contribution_by_{slugify(GENE_BREAKDOWN_GROUP)}.csv")
                display(plot_gene_contribution(gene_contribution, title=f"Relative marker contribution within {GENE_BREAKDOWN_GROUP}"))
                save_current_figure(FIGURES_DIR / f"marker_gene_contribution_by_{slugify(GENE_BREAKDOWN_GROUP)}.png")

                sc.pl.dotplot(
                    adata_sc,
                    var_names=GENE_BREAKDOWN_GENES,
                    groupby=GENE_BREAKDOWN_GROUP,
                    standard_scale="var",
                    dendrogram=False,
                    save=f"_gene_breakdown_{slugify(GENE_BREAKDOWN_GROUP)}.png",
                )
                sc.pl.matrixplot(
                    adata_sc,
                    var_names=GENE_BREAKDOWN_GENES,
                    groupby=GENE_BREAKDOWN_GROUP,
                    standard_scale="var",
                    cmap="mako",
                    save=f"_gene_breakdown_{slugify(GENE_BREAKDOWN_GROUP)}.png",
                )
            else:
                print("Need marker genes plus a valid grouping column such as Renamed_clusternames for cell-type-specific expression views.")
            """
        ),
        code_cell(
            """
            if adata_sc is not None:
                color_columns = [column for column in GROUP_COLUMNS if column in adata_sc.obs.columns]
                color_columns += [column for column in results["score_columns"][:2] if column in adata_sc.obs.columns]
                if color_columns:
                    sc.pl.umap(adata_sc, color=color_columns, ncols=2, save=f"_{slugify(PRIMARY_GROUP or 'group')}_overview.png")
                else:
                    print("Set GROUP_COLUMNS to existing metadata columns before plotting UMAP.")
            """
        ),
        code_cell(
            """
            FEATURE_UMAP_GENES = MARKER_GENES[: min(6, len(MARKER_GENES))]

            if adata_sc is not None and FEATURE_UMAP_GENES:
                sc.pl.umap(adata_sc, color=FEATURE_UMAP_GENES, ncols=3, save=f"_{slugify(PRIMARY_GROUP or 'group')}_marker_genes.png")
            else:
                print("Need marker genes present in the dataset for gene-level UMAP overlays.")
            """
        ),
        code_cell(
            """
            if adata_sc is not None and PRIMARY_GROUP is not None and MARKER_GENES:
                sc.pl.dotplot(adata_sc, var_names=MARKER_GENES, groupby=PRIMARY_GROUP, standard_scale="var", dendrogram=False, save=f"_{slugify(PRIMARY_GROUP)}_markers.png")
                sc.pl.matrixplot(adata_sc, var_names=MARKER_GENES, groupby=PRIMARY_GROUP, standard_scale="var", cmap="mako", save=f"_{slugify(PRIMARY_GROUP)}_markers.png")
            else:
                print("Need valid marker genes and a valid PRIMARY_GROUP for dotplot and matrixplot.")
            """
        ),
        code_cell(
            """
            if adata_sc is not None and PRIMARY_GROUP is not None and results is not None and results["score_columns"]:
                sc.pl.violin(
                    adata_sc,
                    keys=results["score_columns"][: min(3, len(results["score_columns"]))],
                    groupby=PRIMARY_GROUP,
                    stripplot=False,
                    rotation=90,
                    save=f"_{slugify(PRIMARY_GROUP)}_panel_scores.png",
                )
            else:
                print("Need valid score columns and a valid PRIMARY_GROUP for violin plots.")
            """
        ),
        markdown_cell(
            """
            ## Differential Signals
            """
        ),
        code_cell(
            """
            RANK_BY = next((column for column in GROUP_COLUMNS if adata_sc is not None and column in adata_sc.obs.columns), None)

            if adata_sc is not None and RANK_BY is not None:
                rank_df = rank_genes_groups_df(adata_sc, groupby=RANK_BY, repo_root=REPO_ROOT, n_genes=10)
                display(rank_df.head(40))
                save_table(rank_df, RESULTS_DIR / f"rank_genes_by_{slugify(RANK_BY)}.csv", index=False)
                ranked_plot_df = plot_top_ranked_genes(rank_df, top_n=5)
                display(ranked_plot_df)
                save_table(ranked_plot_df, RESULTS_DIR / f"top_ranked_genes_by_{slugify(RANK_BY)}.csv", index=False)
                save_current_figure(FIGURES_DIR / f"top_ranked_genes_by_{slugify(RANK_BY)}.png")
            else:
                print("Set GROUP_COLUMNS so at least one existing obs column is available for marker ranking.")
            """
        ),
        *build_question_specific_cells(question_id),
        markdown_cell(
            """
            ## Interpretation Notes

            This notebook is now scanpy-first. Once you identify useful metadata columns and candidate subpopulations, extend it with:

            - `sc.pl.umap`, `sc.pl.dotplot`, `sc.pl.matrixplot`, or `sc.pl.violin`
            - `sc.tl.rank_genes_groups` on lesion, condition, or cell-state groupings
            - focused filtering to oligodendrocyte-lineage or glial subsets
            - cross-dataset validation using the same panel definitions
            """
        ),
    ]
    return nb


def write_notebook(path: Path, notebook: nbf.NotebookNode) -> None:
    path.write_text(nbf.writes(notebook), encoding="utf-8")


def main() -> None:
    NOTEBOOKS_DIR.mkdir(parents=True, exist_ok=True)
    write_notebook(NOTEBOOKS_DIR / "00_dataset_inventory.ipynb", build_inventory_notebook())
    for question_id in ORDERED_QUESTION_IDS:
        question = QUESTION_BANK[question_id]
        write_notebook(NOTEBOOKS_DIR / f"{question.notebook_slug}.ipynb", build_question_notebook(question_id))


if __name__ == "__main__":
    main()
