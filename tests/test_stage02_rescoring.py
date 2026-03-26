from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from .conftest import run_script, write_text


def test_parse_gnina_outputs_and_integrate_scores(tmp_path: Path) -> None:
    gnina_dir = tmp_path / "gnina_results"
    write_text(
        gnina_dir / "ZINC0001.0.N_gnina.txt",
        """
-----+------------+------------+----------+----------
    1       -5.2      0.6100      6.4100      0.000
        """.strip()
        + "\n",
    )
    write_text(
        gnina_dir / "ZINC0002.0.N_gnina.txt",
        """
-----+------------+------------+----------+----------
    1       -4.8      0.8200      6.8000      0.000
        """.strip()
        + "\n",
    )

    gnina_csv = tmp_path / "outputs" / "gnina_scores.csv"
    parse_result = run_script(
        "02-rescoring/scripts/merge_gnina_results.py",
        [
            "--input-dir",
            str(gnina_dir),
            "--output",
            str(gnina_csv),
        ],
    )
    assert parse_result.returncode == 0, parse_result.stderr

    rf_dir = tmp_path / "rfscore"
    write_text(rf_dir / "ZINC0001.0.N_out.score", "ZINC0001.0.N_out,6.3\n")
    write_text(rf_dir / "ZINC0002.0.N_out.score", "ZINC0002.0.N_out,7.1\n")

    rf_csv = tmp_path / "outputs" / "rf_scores.csv"
    rf_result = run_script(
        "02-rescoring/scripts/merge_rfscore_results.py",
        [
            "--input-dir",
            str(rf_dir),
            "--output",
            str(rf_csv),
        ],
    )
    assert rf_result.returncode == 0, rf_result.stderr

    shortlist_csv = tmp_path / "shortlist.csv"
    with shortlist_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["rank", "zinc_id", "unidock_affinity", "smiles", "subset", "pose_path", "source_pdbqt_path"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "rank": 1,
                "zinc_id": "ZINC0002.0.N",
                "unidock_affinity": -9.1,
                "smiles": "CCN",
                "subset": "zinc-22g",
                "pose_path": "/tmp/ZINC0002.0.N_out.pdbqt",
                "source_pdbqt_path": "archives/a.tgz:subset/ZINC0002.0.N.pdbqt",
            }
        )
        writer.writerow(
            {
                "rank": 2,
                "zinc_id": "ZINC0001.0.N",
                "unidock_affinity": -8.2,
                "smiles": "CCO",
                "subset": "zinc-22g",
                "pose_path": "/tmp/ZINC0001.0.N_out.pdbqt",
                "source_pdbqt_path": "archives/a.tgz:subset/ZINC0001.0.N.pdbqt",
            }
        )

    integrated_csv = tmp_path / "outputs" / "integrated_scores.csv"
    integrate_result = run_script(
        "02-rescoring/scripts/integrate_scores.py",
        [
            "--shortlist-csv",
            str(shortlist_csv),
            "--gnina-csv",
            str(gnina_csv),
            "--rfscore-csv",
            str(rf_csv),
            "--output",
            str(integrated_csv),
        ],
    )
    assert integrate_result.returncode == 0, integrate_result.stderr

    integrated = pd.read_csv(integrated_csv)
    assert list(integrated.columns) == [
        "conformer_id",
        "unidock_affinity",
        "rfscore_v3",
        "gnina_cnnscore",
        "gnina_cnnaffinity",
        "unidock_rank",
        "rfscore_rank",
        "gnina_rank",
        "composite_rank",
    ]
    assert integrated.iloc[0]["conformer_id"] == "ZINC0002.0.N"
    assert integrated["composite_rank"].min() == integrated.iloc[0]["composite_rank"]
