from pathlib import Path
from typing import List, Optional, Tuple

from rdkit import Chem

_NAME_TO_SMILES = None


def _load_name_map() -> dict:
    """Load name→SMILES from local CSV (reverse of smiles→name in common_names.csv)."""
    mapping = {}
    csv_path = Path(__file__).resolve().parent / "name_map.csv"
    if not csv_path.exists():
        return mapping
    import csv
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header
        for row in reader:
            if len(row) >= 2:
                smiles, name = row[0].strip(), row[1].strip()
                # Normalize name to lower-case for case-insensitive lookup
                mapping[name.lower()] = smiles
                # Also index without common punctuation
                simple = name.lower().replace("-", " ").replace(",", "")
                if simple != name.lower():
                    mapping[simple] = smiles
    return mapping


def canonicalize_smiles(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))


def is_smiles(candidate: str) -> bool:
    return canonicalize_smiles(candidate) is not None


def _resolve_via_pubchem(name: str) -> Optional[str]:
    """Try to resolve a molecule name to SMILES using cirpy or pubchempy."""
    try:
        import cirpy
        smiles = cirpy.resolve(name, "smiles")
        if smiles and canonicalize_smiles(smiles):
            return smiles
    except Exception:
        pass
    try:
        from pubchempy import get_compounds
        compounds = get_compounds(name, "name", timeout=5)
        if compounds:
            return compounds[0].canonical_smiles
    except Exception:
        pass
    return None


def resolve_identifier(identifier: str) -> Tuple[Optional[str], str]:
    canonical = canonicalize_smiles(identifier)
    if canonical is not None:
        return canonical, "smiles"

    # Try local name→SMILES map first
    global _NAME_TO_SMILES
    if _NAME_TO_SMILES is None:
        _NAME_TO_SMILES = _load_name_map()
    smiles = _NAME_TO_SMILES.get(identifier.lower())
    if smiles and canonicalize_smiles(smiles):
        return smiles, "local"

    # Fallback: online resolution
    smiles = _resolve_via_pubchem(identifier)
    if smiles:
        canonical = canonicalize_smiles(smiles)
        if canonical:
            return canonical, "pubchem"

    return None, "unresolved"


def resolve_identifier_list(identifiers: List[str]) -> Tuple[List[str], List[str]]:
    resolved = []
    unresolved = []
    for identifier in identifiers:
        smiles, _ = resolve_identifier(identifier)
        if smiles is None:
            unresolved.append(identifier)
        else:
            resolved.append(smiles)
    return resolved, unresolved
