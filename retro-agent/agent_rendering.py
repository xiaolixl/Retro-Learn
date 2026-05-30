from pathlib import Path
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D


def _draw_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for rendering: {smiles}")
    AllChem.Compute2DCoords(mol)
    return mol


def _save_png(mol: Chem.Mol, output_path: Path, size=(420, 280)) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image = Draw.MolToImage(mol, size=size)
    image.save(output_path)
    return str(output_path)


def _save_svg(mol: Chem.Mol, output_path: Path, size=(300, 200)) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    output_path.write_text(drawer.GetDrawingText(), encoding="utf-8")
    return str(output_path)


def _save_molecule(smiles: str, output_dir: Path, name: str) -> Dict[str, str]:
    mol = _draw_mol(smiles)
    return {
        "png": _save_png(mol, output_dir / f"{name}.png"),
        "svg": _save_svg(mol, output_dir / f"{name}.svg"),
    }


def _save_grid_svg(
    mols: List[Chem.Mol],
    legends: List[str],
    output_path: Path,
    mol_size=(220, 160),
) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cols = min(3, max(1, len(mols)))
    rows = (len(mols) + cols - 1) // cols
    total_w = cols * mol_size[0]
    total_h = rows * mol_size[1]
    drawer = rdMolDraw2D.MolDraw2DSVG(total_w, total_h, mol_size[0], mol_size[1])
    drawer.DrawMolecules(mols, legends=legends)
    drawer.FinishDrawing()
    output_path.write_text(drawer.GetDrawingText(), encoding="utf-8")
    return str(output_path)


def _save_grid_png(
    mols: List[Chem.Mol],
    legends: List[str],
    output_path: Path,
    mol_size=(220, 160),
) -> str:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    image = Draw.MolsToGridImage(
        mols,
        molsPerRow=min(3, max(1, len(mols))),
        subImgSize=mol_size,
        legends=legends,
        useSVG=False,
    )
    image.save(output_path)
    return str(output_path)


def render_agent_artifacts(result: Dict[str, Any], output_dir: Path) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = result["data"]
    target_smiles = data["target_molecule"]["smiles"]
    target = _save_molecule(target_smiles, output_dir, "target")

    if data["mode"] == "single_step":
        route_images = []
        routes = data.get("recommended_routes") or data.get("retrosynthesis_routes", [])
        for route in routes:
            route_dir = output_dir / f"route_{route['route_rank']}"
            reactants = []
            reactant_mols = []
            legends = []
            for idx, reactant in enumerate(route["reactants"], start=1):
                imgs = _save_molecule(reactant["smiles"], route_dir, f"reactant_{idx}")
                reactants.append({"smiles": reactant["smiles"], **imgs})
                reactant_mols.append(_draw_mol(reactant["smiles"]))
                legends.append(f"Reactant {idx}")
            route_grid_png = _save_grid_png(reactant_mols, legends, route_dir / "route_grid.png")
            route_grid_svg = _save_grid_svg(reactant_mols, legends, route_dir / "route_grid.svg")
            route_images.append({
                "route_rank": route["route_rank"],
                "route_grid_png": route_grid_png,
                "route_grid_svg": route_grid_svg,
                "reactants": reactants,
            })
        return {"target": target, "routes": route_images}

    # multi-step
    recommended_route = data["recommended_route"]
    if recommended_route is None:
        return {"target": target, "routes": []}

    step_images = []
    for step in recommended_route["steps"]:
        step_dir = output_dir / f"step_{step['step_number']}"
        expanded = _save_molecule(step["expanded_smiles"], step_dir, "expanded_target")
        # Parse reactants from expanded_smiles (may be "A.B" for multi-reactant steps)
        reactants = []
        reactant_mols = []
        legends = []
        expanded_parts = step["expanded_smiles"].split(".")
        for idx, r_smiles in enumerate(expanded_parts, start=1):
            imgs = _save_molecule(r_smiles, step_dir, f"reactant_{idx}")
            reactants.append({"smiles": r_smiles, **imgs})
            reactant_mols.append(_draw_mol(r_smiles))
            legends.append(f"Step {step['step_number']} R{idx}")
        grid_png = _save_grid_png(reactant_mols, legends, step_dir / "step_grid.png")
        grid_svg = _save_grid_svg(reactant_mols, legends, step_dir / "step_grid.svg")
        step_images.append({
            "step_number": step["step_number"],
            "expanded_smiles": step["expanded_smiles"],
            "expanded": expanded,
            "step_grid_png": grid_png,
            "step_grid_svg": grid_svg,
            "reactants": reactants,
        })

    final_leaf_dir = output_dir / "final_leaf_reactants"
    final_leaf_images = []
    leaf_reactants = recommended_route.get("leaf_reactants") or recommended_route.get("final_leaf_reactants", [])
    for idx, reactant in enumerate(leaf_reactants, start=1):
        imgs = _save_molecule(reactant["smiles"], final_leaf_dir, f"leaf_{idx}")
        final_leaf_images.append({"smiles": reactant["smiles"], **imgs})

    return {
        "target": target,
        "routes": step_images,
        "final_leaf_reactants": final_leaf_images,
    }
