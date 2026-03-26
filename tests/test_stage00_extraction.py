from __future__ import annotations

import csv
import tarfile
from pathlib import Path

from .conftest import run_script


def _make_archive(archive_path: Path, members: dict[str, str]) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "w:gz") as tar:
        for internal_path, content in members.items():
            tmp = archive_path.parent / Path(internal_path).name
            tmp.write_text(content)
            tar.add(tmp, arcname=internal_path)
            tmp.unlink()


def test_extract_pdbqt_smoke(tmp_path: Path) -> None:
    data_root = tmp_path / "data_root"
    archive = data_root / "subset" / "batch.tgz"
    _make_archive(
        archive,
        {
            "subset/ZINC0001.0.N.pdbqt": "REMARK Name = ZINC0001.0.N\nATOM\nENDMDL\n",
            "subset/ZINC0002.0.N.pdbqt": "REMARK Name = ZINC0002.0.N\nATOM\nENDMDL\n",
        },
    )

    input_csv = tmp_path / "inputs.csv"
    with input_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["zinc_id", "smiles", "pdbqt_path"])
        writer.writeheader()
        writer.writerow(
            {
                "zinc_id": "ZINC0001.0.N",
                "smiles": "CCO",
                "pdbqt_path": "subset/batch.tgz:subset/ZINC0001.0.N.pdbqt",
            }
        )
        writer.writerow(
            {
                "zinc_id": "ZINC0002.0.N",
                "smiles": "CCN",
                "pdbqt_path": "subset/batch.tgz:subset/ZINC0002.0.N.pdbqt",
            }
        )

    output_dir = tmp_path / "ligands"
    ligand_list = tmp_path / "outputs" / "ligand_files.txt"
    result = run_script(
        "00-data-extraction/scripts/extract_pdbqt.py",
        [
            "--input",
            str(input_csv),
            "--output_dir",
            str(output_dir),
            "--data_root",
            str(data_root),
            "--ligand_list",
            str(ligand_list),
            "--ncpus",
            "1",
        ],
    )

    assert result.returncode == 0, result.stderr
    extracted = sorted(p.name for p in output_dir.glob("*.pdbqt"))
    assert extracted == ["ZINC0001.0.N.pdbqt", "ZINC0002.0.N.pdbqt"]
    assert ligand_list.read_text().count(".pdbqt") == 2
