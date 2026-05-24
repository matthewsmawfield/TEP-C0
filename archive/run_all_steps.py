#!/usr/bin/env python3
"""TEP-C0 master pipeline driver.

This wrapper delegates to the canonical pipeline in ``scripts/run_pipeline.py``.
It keeps older commands working without maintaining a second, divergent step
list.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_pipeline import run_pipeline


def main() -> int:
    results = run_pipeline()
    failed = [name for name, result in results.items() if "error" in result]
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
