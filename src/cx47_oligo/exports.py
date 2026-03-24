from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


def slugify(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "value"


def dataset_slug(dataset_path: str | Path | None) -> str:
    if dataset_path is None:
        return "no_dataset"
    return slugify(Path(dataset_path).stem)


def ensure_results_dir(repo_root: str | Path, notebook_slug: str, dataset_path: str | Path | None = None) -> Path:
    root = Path(repo_root)
    output_dir = root / "results" / slugify(notebook_slug) / dataset_slug(dataset_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def save_table(frame: pd.DataFrame, path: str | Path, index: bool = True) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, index=index)
    return output_path


def save_text(text: str, path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    return output_path


def save_json(data: dict[str, Any], path: str | Path) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    return output_path


def save_current_figure(path: str | Path, close: bool = True, dpi: int = 240) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.gcf().savefig(output_path, dpi=dpi, bbox_inches="tight")
    if close:
        plt.close(plt.gcf())
    return output_path
