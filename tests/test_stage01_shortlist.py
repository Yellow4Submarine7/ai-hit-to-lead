from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from .conftest import run_script, write_text


def test_export_results_builds_shortlist_contract(tmp_path: Path) -> None:
    results_dir = tmp_path / "results"
    write_text(
        results_dir / "ZINC0001.0.N_out.pdbqt",
        "REMARK VINA RESULT:    -8.20      0.000      0.000\nREMARK Name = ZINC0001.0.N\n",
    )
    write_text(
        results_dir / "ZINC0002.0.N_out.pdbqt",
        "REMARK VINA RESULT:    -9.10      0.000      0.000\nREMARK Name = ZINC0002.0.N\n",
    )

    input_csv = tmp_path / "smiles.csv"
    with input_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["zinc_id", "smiles", "pdbqt_path", "subset"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "zinc_id": "ZINC0001.0.N",
                "smiles": "CCO",
                "pdbqt_path": "archives/a.tgz:subset/ZINC0001.0.N.pdbqt",
                "subset": "zinc-22g",
            }
        )
        writer.writerow(
            {
                "zinc_id": "ZINC0002.0.N",
                "smiles": "CCN",
                "pdbqt_path": "archives/a.tgz:subset/ZINC0002.0.N.pdbqt",
                "subset": "zinc-22g",
            }
        )

    output_csv = tmp_path / "outputs" / "shortlist.csv"
    manifest_output = tmp_path / "outputs" / "shortlist_manifest.txt"
    result = run_script(
        "01-gpu-docking/scripts/export_results.py",
        [
            "--results-dir",
            str(results_dir),
            "--smiles-csv",
            str(input_csv),
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--manifest-output",
            str(manifest_output),
            "--top-n",
            "1",
            "--ncpus",
            "1",
            "--batch-size",
            "10",
        ],
    )

    assert result.returncode == 0, result.stderr
    exported = pd.read_csv(output_csv)
    assert list(exported.columns) == [
        "rank",
        "zinc_id",
        "unidock_affinity",
        "smiles",
        "subset",
        "pose_path",
        "source_pdbqt_path",
    ]
    assert exported.iloc[0]["zinc_id"] == "ZINC0002.0.N"
    assert exported.iloc[0]["rank"] == 1
    assert manifest_output.read_text().strip().endswith("ZINC0002.0.N_out.pdbqt")
