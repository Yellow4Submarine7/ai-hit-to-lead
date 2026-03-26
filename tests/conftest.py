from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON = REPO_ROOT / ".venv" / "bin" / "python"


def run_script(relative_path: str, args: Iterable[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    script_path = REPO_ROOT / relative_path
    return subprocess.run(
        [str(PYTHON), str(script_path), *args],
        cwd=str(cwd or REPO_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )


def write_text(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path
