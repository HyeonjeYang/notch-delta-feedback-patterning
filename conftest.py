"""Ensures the repository root is importable (model, geometry, stability,
control, simulation, figures packages) regardless of pytest's invocation
directory."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
