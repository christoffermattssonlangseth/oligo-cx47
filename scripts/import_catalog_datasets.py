from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from cx47_oligo.h5ad_tools import dataset_catalog
from cx47_oligo.source_import import import_catalog_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Import configured matrix-plus-metadata datasets into local h5ad files.")
    parser.add_argument(
        "--dataset-id",
        action="append",
        dest="dataset_ids",
        help="Limit import to one or more dataset IDs from config/dataset_catalog.csv.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing imported h5ad files.",
    )
    return parser.parse_args()


def select_rows(catalog: pd.DataFrame, dataset_ids: list[str] | None) -> pd.DataFrame:
    rows = catalog.dropna(subset=["local_path", "source_expression_path", "source_metadata_path"]).copy()
    if dataset_ids:
        rows = rows[rows["dataset_id"].isin(dataset_ids)].copy()
    return rows


def main() -> None:
    args = parse_args()
    catalog = dataset_catalog(REPO_ROOT / "config" / "dataset_catalog.csv")
    rows = select_rows(catalog, args.dataset_ids)
    if rows.empty:
        raise SystemExit("No matching importable datasets were found in config/dataset_catalog.csv.")

    for _, row in rows.iterrows():
        output_path = REPO_ROOT / row["local_path"]
        if output_path.exists() and not args.force:
            print(f"Skipping {row['dataset_id']} because {output_path} already exists. Use --force to rebuild.")
            continue

        print(f"Importing {row['dataset_id']} -> {output_path}")
        written_path = import_catalog_dataset(row, REPO_ROOT)
        print(f"Wrote {written_path}")


if __name__ == "__main__":
    main()
