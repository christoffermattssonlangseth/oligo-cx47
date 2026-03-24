from __future__ import annotations

import os
from pathlib import Path


def configure_scanpy_environment(repo_root: str | Path | None = None) -> dict[str, Path]:
    root = Path(repo_root).resolve() if repo_root is not None else Path.cwd().resolve()
    cache_root = root / ".cache" / "scanpy"
    numba_cache = cache_root / "numba"
    mpl_cache = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"

    for directory in (numba_cache, mpl_cache, xdg_cache):
        directory.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("NUMBA_CACHE_DIR", str(numba_cache))
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))
    os.environ.setdefault("LOKY_MAX_CPU_COUNT", str(max(1, os.cpu_count() or 1)))
    return {
        "cache_root": cache_root,
        "numba_cache": numba_cache,
        "mpl_cache": mpl_cache,
        "xdg_cache": xdg_cache,
    }


def import_scanpy(repo_root: str | Path | None = None):
    cache_paths = configure_scanpy_environment(repo_root=repo_root)
    import numba
    from numba.core import config as numba_config

    numba_config.CACHE_DIR = str(cache_paths["numba_cache"])
    import scanpy as sc

    sc.settings.verbosity = 2
    sc.settings.n_jobs = 1
    return sc
