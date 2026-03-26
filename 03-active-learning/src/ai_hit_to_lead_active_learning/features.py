"""
Feature extraction helpers for the public Teacher ensemble.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from sklearn.preprocessing import StandardScaler


PHYSICHEM_DESCRIPTORS = [
    ("MolWt", Descriptors.MolWt),
    ("LogP", Descriptors.MolLogP),
    ("TPSA", Descriptors.TPSA),
    ("NumHDonors", Descriptors.NumHDonors),
    ("NumHAcceptors", Descriptors.NumHAcceptors),
    ("NumRotatableBonds", Descriptors.NumRotatableBonds),
    ("HeavyAtomCount", Descriptors.HeavyAtomCount),
    ("FractionCSP3", Descriptors.FractionCSP3),
]


def smiles_to_mol(smiles: str) -> Optional[Chem.Mol]:
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def compute_morgan_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048) -> np.ndarray:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    arr = np.zeros(n_bits, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def compute_physichem_descriptors(mol: Chem.Mol) -> np.ndarray:
    values = []
    for _, func in PHYSICHEM_DESCRIPTORS:
        try:
            value = func(mol)
            if value is None or np.isnan(value) or np.isinf(value):
                value = 0.0
        except Exception:
            value = 0.0
        values.append(value)
    return np.array(values, dtype=np.float32)


def extract_features(
    smiles_list: List[str],
    fit_scaler: bool = False,
    scaler: Optional[StandardScaler] = None,
    radius: int = 2,
    n_bits: int = 2048,
) -> Tuple[np.ndarray, Optional[StandardScaler]]:
    n_items = len(smiles_list)
    morgan = np.zeros((n_items, n_bits), dtype=np.float32)
    physichem = np.zeros((n_items, len(PHYSICHEM_DESCRIPTORS)), dtype=np.float32)

    for index, smiles in enumerate(smiles_list):
        mol = smiles_to_mol(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES at index {index}: {smiles}")
        morgan[index] = compute_morgan_fingerprint(mol, radius=radius, n_bits=n_bits)
        physichem[index] = compute_physichem_descriptors(mol)

    if fit_scaler:
        scaler = StandardScaler()
        physichem = scaler.fit_transform(physichem)
    elif scaler is not None:
        physichem = scaler.transform(physichem)

    return np.concatenate([morgan, physichem], axis=1), scaler
