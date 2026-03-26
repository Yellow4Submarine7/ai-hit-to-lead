from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from .conftest import run_script


def _write_experimental_csv(path: Path) -> None:
    rows = [
        ("MCULE-1", "CCO", -1.2, "scaffold_a", "/tmp/MCULE-1.pdbqt"),
        ("MCULE-2", "CCN", -0.8, "scaffold_a", "/tmp/MCULE-2.pdbqt"),
        ("MCULE-3", "CCC", -0.2, "scaffold_b", "/tmp/MCULE-3.pdbqt"),
        ("MCULE-4", "CCCl", 0.1, "scaffold_b", "/tmp/MCULE-4.pdbqt"),
        ("MCULE-5", "CCBr", 0.6, "scaffold_c", "/tmp/MCULE-5.pdbqt"),
        ("MCULE-6", "CCF", 1.0, "scaffold_c", "/tmp/MCULE-6.pdbqt"),
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mcule_id", "smiles_canonical", "binding_level_nm", "scaffold", "pdbqt_path"])
        writer.writerows(rows)


def _write_candidates_csv(path: Path) -> None:
    rows = [
        ("ZINC0001.0.N", "CCO", "/tmp/ZINC0001.0.N_out.pdbqt"),
        ("ZINC0002.0.N", "CCN", "/tmp/ZINC0002.0.N_out.pdbqt"),
        ("ZINC0003.0.N", "CCC", "/tmp/ZINC0003.0.N_out.pdbqt"),
        ("ZINC0004.0.N", "CCCl", "/tmp/ZINC0004.0.N_out.pdbqt"),
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand_id", "smiles_canonical", "pdbqt_path"])
        writer.writerows(rows)


def test_train_teacher_and_generate_pseudo_labels(tmp_path: Path) -> None:
    experimental_csv = tmp_path / "experimental.csv"
    candidates_csv = tmp_path / "candidates.csv"
    _write_experimental_csv(experimental_csv)
    _write_candidates_csv(candidates_csv)

    checkpoint = tmp_path / "checkpoints" / "teacher.pkl"
    metrics_csv = tmp_path / "results" / "teacher_metrics.csv"
    train_result = run_script(
        "03-active-learning/scripts/train_teacher.py",
        [
            "--experimental-csv",
            str(experimental_csv),
            "--checkpoint-output",
            str(checkpoint),
            "--metrics-output",
            str(metrics_csv),
            "--seeds",
            "11,12",
        ],
    )
    assert train_result.returncode == 0, train_result.stderr
    assert checkpoint.exists()
    metrics = pd.read_csv(metrics_csv)
    assert {"seed", "spearman", "train_size", "val_size"} <= set(metrics.columns)

    pseudo_dir = tmp_path / "pseudo"
    pseudo_result = run_script(
        "03-active-learning/scripts/generate_pseudo_labels.py",
        [
            "--teacher-checkpoint",
            str(checkpoint),
            "--experimental-csv",
            str(experimental_csv),
            "--candidates-csv",
            str(candidates_csv),
            "--receptor-path",
            "receptor/wwp2_c2.pdbqt",
            "--output-dir",
            str(pseudo_dir),
            "--k-values",
            "2,4",
        ],
    )
    assert pseudo_result.returncode == 0, pseudo_result.stderr
    assert (pseudo_dir / "teacher_predictions.csv").exists()
    assert (pseudo_dir / "train_K2.types").exists()
    assert (pseudo_dir / "train_K4.types").exists()
    assert (pseudo_dir / "val.types").exists()


def test_score_candidates_with_precomputed_predictions(tmp_path: Path) -> None:
    candidates_csv = tmp_path / "candidates.csv"
    _write_candidates_csv(candidates_csv)

    predictions_csv = tmp_path / "predictions.csv"
    with predictions_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand_id", "affinity_pred"])
        writer.writerow(["ZINC0001.0.N", -6.1])
        writer.writerow(["ZINC0002.0.N", -7.3])
        writer.writerow(["ZINC0003.0.N", -5.2])
        writer.writerow(["ZINC0004.0.N", -8.0])

    ranked_csv = tmp_path / "ranked.csv"
    result = run_script(
        "03-active-learning/scripts/score_candidates.py",
        [
            "--candidates-csv",
            str(candidates_csv),
            "--predictions-csv",
            str(predictions_csv),
            "--output",
            str(ranked_csv),
        ],
    )
    assert result.returncode == 0, result.stderr
    ranked = pd.read_csv(ranked_csv)
    assert list(ranked.columns) == ["rank", "ligand_id", "smiles_canonical", "student_score"]
    assert ranked.iloc[0]["ligand_id"] == "ZINC0004.0.N"
