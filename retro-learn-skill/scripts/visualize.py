#!/usr/bin/env python3
"""
SimpRetro result visualizer.
Converts retro_result.json into an HTML page with molecular structure drawings.
Forward-synthesis direction · single-arrow layout · 60% SVG size.

Usage:
  python visualize.py retro_result.json                    # output: retro_result_view.html
  python visualize.py retro_result.json -o custom.html
"""

import argparse
import csv
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def smiles_to_svg_full(smiles, width=240, height=180):
    """Convert SMILES to inline SVG string, keeping ALL elements intact."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f'<div style="color:#999;font-size:11px;">Invalid: {smiles}</div>'
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Remove XML PI only — keep everything else (atom labels, double bonds, backgrounds)
    svg = svg.replace('<?xml version="1.0"?>', '')
    svg = svg.replace('<svg', '<svg class="mol-svg"')
    return svg


def esc(text):
    """Escape HTML special characters."""
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _load_common_names():
    """Load SMILES→common_name mapping from CSV file."""
    names = {}
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "common_names.csv")
    if not os.path.exists(csv_path):
        return names
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader, None)  # skip header
        for row in reader:
            if len(row) >= 2:
                names[row[0].strip()] = row[1].strip()
    return names


_COMMON_NAMES = _load_common_names()


def _common_name(smiles):
    """Return common name for a SMILES, or empty string if unknown."""
    if not smiles:
        return ""
    # Direct lookup
    if smiles in _COMMON_NAMES:
        return _COMMON_NAMES[smiles]
    # Try canonicalizing via RDKit
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            cansmi = Chem.MolToSmiles(mol)
            if cansmi in _COMMON_NAMES:
                return _COMMON_NAMES[cansmi]
    except Exception:
        pass
    return ""


def classify_reaction(template_str, condition_list, reaction_type=None):
    """Return reaction type label.

    Prefers pre-classified reaction_type from template cache (injected by
    run_retro.py). Falls back to empty string if not classified.
    """
    if reaction_type:
        return reaction_type
    return ""


def _rxn_type_class(label):
    ll = label.lower()
    if "diels-alder" in ll:
        return "rxn-type-da"
    if "esterification" in ll or "amidation" in ll or "acylation" in ll:
        return "rxn-type-ester"
    if "coupling" in ll or "cross-coupling" in ll:
        return "rxn-type-coupling"
    if "elimination" in ll or "dehydrohalogenation" in ll:
        return "rxn-type-elim"
    if "reduct" in ll or "hydrogenation" in ll:
        return "rxn-type-red"
    if "oxid" in ll:
        return "rxn-type-ox"
    if "halogen" in ll or "bromin" in ll or "chlorin" in ll or "iodin" in ll:
        return "rxn-type-hal"
    if "hydrolysis" in ll:
        return "rxn-type-hydro"
    if "substitution" in ll or "sn1" in ll or "sn2" in ll:
        return "rxn-type-sub"
    if "addition" in ll or "grignard" in ll or "wittig" in ll or "aldol" in ll:
        return "rxn-type-add"
    if "alkylation" in ll or "alkylat" in ll:
        return "rxn-type-alk"
    return "rxn-type-default"


def _dedup_conditions(cond_list):
    """Deduplicate reaction conditions, preserving order."""
    if not cond_list:
        return ""
    seen = set()
    unique = []
    for c in cond_list:
        c = str(c).strip()
        if c and c not in seen:
            seen.add(c)
            unique.append(c)
    return ", ".join(unique)


def _make_arrow_section(conditions, rxn_label, rxn_cls):
    """Build the single-arrow HTML with conditions above, reaction type below."""
    cond_html = esc(conditions) if conditions else ""
    return f'''<div class="arrow-section">
        <span class="cond">{cond_html}</span>
        <span class="rxn-type {rxn_cls}">{esc(rxn_label)}</span>
        <span class="arrow-line"></span>
    </div>'''


def _make_mol_box(svg, smiles, name="", stock=None):
    """Build a molecule display box, with common name if known."""
    stock_tag = ""
    if stock is not None:
        if stock:
            stock_tag = '<span class="stock-in">In Stock</span>'
    name_html = f'<span class="hl">{esc(name)}</span><br>' if name else ""
    cname = _common_name(smiles)
    cname_html = f'<br><span class="cname">{esc(cname)}</span>' if cname else ""
    return f'''<div class="mol-wrap">
        {svg}
        <div class="ml">{name_html}{esc(smiles)}{cname_html} {stock_tag}</div>
    </div>'''


def _make_target_mol_box(svg, smiles, name="", mw=""):
    """Build a target molecule display box with distinctive 'target' label.

    Same layout as regular mol-box (name → SMILES) but with indigo border,
    indigo background on SVG, and a 'target' pill badge.
    """
    name_html = f'<span class="hl">{esc(name)}</span><br>' if name else ""
    mw_html = f'<br>MW: {esc(mw)}' if mw else ""
    return f'''<div class="mol-wrap mol-target">
        {svg}
        <div class="ml">{name_html}{esc(smiles)}{mw_html}
            <span class="target-label">target</span>
        </div>
    </div>'''


def generate_flowchart_html(data, output_path):
    """Generate a flowchart-style HTML visualization with forward-direction single-arrow layout."""
    target = data["data"]["target_molecule"]
    mode = data["data"].get("mode", "single_step")
    target_svg = smiles_to_svg_full(target["smiles"], 200, 160)
    target_name = _common_name(target["smiles"]) or "Target"
    target_mw = target.get("molecular_weight", "")

    html_css = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SimpRetro Retrosynthesis Result</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, 'Segoe UI', sans-serif; background: #f8f8f6; color: #2c2c2a; padding: 24px; max-width: 850px; margin: 0 auto; }
  h1 { font-size: 18px; font-weight: 500; margin-bottom: 20px; color: #26215C; }
  .target-card { background: #EEEDFE; border: 1.5px solid #534AB7; border-radius: 12px; padding: 14px; text-align: center; margin-bottom: 20px; }
  .target-card .label { font-size: 13px; font-weight: 500; color: #26215C; margin-bottom: 2px; }
  .target-card .smiles { font-family: 'Courier New', monospace; font-size: 12px; color: #534AB7; }
  .target-card .mw { font-size: 11px; color: #3C3489; }
  .target-card .mol-container { display: flex; justify-content: center; margin: 6px 0; }
  .target-card .mol-svg { max-width: 140px; height: auto; }

  .route-card { background: #fff; border-radius: 10px; padding: 16px 20px; margin-bottom: 16px; border: 1.5px solid; }
  .route-best { background: #e1f5ee; border-color: #0F6E56; }
  .route-normal { background: #f1efe8; border-color: #b4b2a9; }
  .route-header { display: flex; align-items: center; gap: 8px; margin-bottom: 10px; flex-wrap: wrap; }
  .route-badge { display: inline-block; padding: 2px 8px; border-radius: 5px; font-size: 11px; font-weight: 500; color: #fff; }
  .route-best .route-badge { background: #0F6E56; }
  .route-normal .route-badge { background: #5F5E5A; }
  .route-score { font-size: 12px; font-weight: 500; }
  .route-best .route-score { color: #04342C; }
  .route-normal .route-score { color: #444441; }
  .best-tag { display: inline-block; padding: 2px 8px; border-radius: 8px; background: rgba(29,158,117,0.15); font-size: 10px; font-weight: 500; color: #0F6E56; }

  /* Single-arrow reaction step */
  .step { display: flex; align-items: center; justify-content: center; gap: 0; padding: 6px 0; flex-wrap: wrap; }
  .reactants { display: flex; align-items: center; gap: 6px; flex-wrap: wrap; justify-content: center; }
  .mol-wrap { display: flex; flex-direction: column; align-items: center; }
  .mol-wrap .mol-svg { max-width: 100px; height: auto; }
  .mol-wrap .ml { font-size: 10px; color: #5F5E5A; text-align: center; margin-top: 2px; line-height: 1.3; }
  .mol-wrap .ml .hl { font-weight: 500; color: #26215C; }
  .mol-wrap .ml .stock-in, .stock-in { display: inline-block; padding: 1px 5px; border-radius: 3px; font-size: 9px; background: #e1f5ee; color: #085041; }
  .mol-wrap .ml .stock-out, .stock-out { display: inline-block; padding: 1px 5px; border-radius: 3px; font-size: 9px; background: #fcebeb; color: #791f1f; }
  .mol-wrap .ml .cname, .cname { font-size: 9px; color: #6b6b6b; font-style: italic; }
  .plus { font-size: 14px; color: #888; margin: 0 2px; }

  /* Arrow with conditions above */
  .arrow-section { display: flex; flex-direction: column; align-items: center; min-width: 100px; flex-shrink: 0; }
  .arrow-section .cond { font-size: 11px; font-weight: 500; color: #3C3489; text-align: center; line-height: 1.3; }
  .arrow-section .rxn-type { font-size: 9px; margin-top: 2px; padding: 1px 5px; border-radius: 3px; font-weight: 500; display: inline-block; }
  .arrow-section .arrow-line { display: block; width: 100%; height: 2px; background: #534AB7; margin: 4px 0; position: relative; }
  .arrow-section .arrow-line::after { content: ''; position: absolute; right: -2px; top: -5px; border-left: 8px solid #534AB7; border-top: 6px solid transparent; border-bottom: 6px solid transparent; }

  .rxn-type-elim { background: #e8f4fd; color: #1a5c8a; }
  .rxn-type-red { background: #e8fde8; color: #1a6b1a; }
  .rxn-type-ox { background: #fde8e8; color: #8a1a1a; }
  .rxn-type-hal { background: #f3e8fd; color: #5c1a8a; }
  .rxn-type-da { background: #fffbe0; color: #6b5a00; }
  .rxn-type-ester { background: #fde8f3; color: #8a1a5c; }
  .rxn-type-coupling { background: #fff0e0; color: #8a4a1a; }
  .rxn-type-hydro { background: #e0f8f8; color: #1a6b6b; }
  .rxn-type-sub { background: #e0f0f8; color: #1a4a6b; }
  .rxn-type-add { background: #e8fde4; color: #2a6b1a; }
  .rxn-type-alk { background: #f0e8fd; color: #4a1a8a; }
  .rxn-type-default { background: #f0f0f0; color: #555; }

  .step-chain { display: block; padding: 4px 0; }
  .step-row { display: block; margin-bottom: 4px; }
  .chain-mol { display: flex; flex-direction: column; align-items: center; }
  .chain-mol .mol-svg { max-width: 90px; height: auto; }
  .chain-mol .cl { font-size: 9px; color: #5F5E5A; text-align: center; margin-top: 1px; line-height: 1.2; }
  .chain-mol .cl .cname { font-size: 8px; color: #6b6b6b; font-style: italic; }

  .step-label { font-size: 10px; color: #aaa; margin: 10px 0 2px 0; padding-left: 6px; letter-spacing: 0.3px; }
  .step-label:first-child { margin-top: 4px; }
  .info-line { font-size: 11px; color: #555; margin-top: 6px; line-height: 1.6; }
  .footer { font-size: 11px; color: #999; text-align: center; margin-top: 20px; padding: 12px; border-top: 1px solid #eee; }

  /* Target molecule distinctive styling */
  .mol-target { position: relative; }
  .mol-target .mol-svg { border: 2px solid #4338CA; border-radius: 6px; background: #EEEDFE; }
  .target-label { display: inline-block; margin-left: 6px; padding: 2px 8px; border-radius: 4px; font-size: 10px; font-weight: 600; color: #fff; background: #4338CA; letter-spacing: 0.5px; vertical-align: middle; }
</style>
</head>
<body>
<h1>Synthetic Route Analysis</h1>
"""

    parts = [html_css]

    # Target card — uses _make_target_mol_box for name + SMILES + target label
    parts.append(f"""<div class="target-card">
  {_make_target_mol_box(target_svg, target['smiles'], target_name, str(target_mw))}
</div>""")

    if mode == "single_step":
        routes = data["data"].get("retrosynthesis_routes", [])
        for route in routes:
            rank = route["route_rank"]
            score = route["score"]
            is_best = rank == 1 and score > 0
            card_cls = "route-best" if is_best else "route-normal"
            cond = _dedup_conditions(route.get("reaction_condition", []))
            rxn_label = classify_reaction(route.get("reaction_template", ""), route.get("reaction_condition", []), route.get("reaction_type"))
            rxn_cls = _rxn_type_class(rxn_label)

            parts.append(f'<div class="route-card {card_cls}">')
            parts.append(f'  <div class="route-header">')
            parts.append(f'    <span class="route-badge">Route {rank}</span>')
            parts.append(f'    <span class="route-score">Score: {score:.4f}</span>')
            if is_best:
                parts.append(f'    <span class="best-tag">Best</span>')
            parts.append(f'  </div>')
            parts.append(f'  <div class="step">')

            # Reactants on the left (forward direction)
            parts.append(f'    <div class="reactants">')
            for ri, reactant in enumerate(route.get("reactants", [])):
                r_svg = smiles_to_svg_full(reactant["smiles"], 200, 160)
                parts.append(f'      {_make_mol_box(r_svg, reactant["smiles"], "", stock=reactant.get("in_stock"))}')
                if ri < len(route["reactants"]) - 1:
                    parts.append(f'      <span class="plus">+</span>')
            parts.append(f'    </div>')

            # Single arrow with conditions on top
            parts.append(f'    {_make_arrow_section(cond, rxn_label, rxn_cls)}')

            # Target on the right — with target label
            parts.append(f'    {_make_target_mol_box(target_svg, target["smiles"], target_name, str(target_mw))}')

            parts.append(f'  </div>')
            parts.append(f'</div>')

    else:
        # Multi-step — show up to 3 routes from all_routes
        all_routes = data["data"].get("all_routes", [])
        recommended = data["data"].get("recommended_route")
        if not all_routes and recommended:
            # Build synthetic entry from recommended_route
            all_routes = [{
                "route_rank": 1,
                "route_score": recommended.get("route_score", 0),
                "steps": recommended.get("actual_steps", 0),
                "leaf_reactants": recommended.get("leaf_reactants", []),
                "steps_history": recommended.get("steps", []),
            }]

        shown = 0
        for route in all_routes[:3]:
            steps_history = route.get("steps_history", [])
            rank = route.get("route_rank", shown + 1)
            score = route.get("route_score", 0)
            step_count = route.get("steps", len(steps_history))
            is_best = rank == 1 and score > 0
            card_cls = "route-best" if is_best else "route-normal"

            chain_parts = []
            # Build stock lookup from leaf_reactants (canonical SMILES → in_stock)
            leaves = route.get("leaf_reactants", [])
            stock_lookup = {}
            for leaf in leaves:
                try:
                    mol = Chem.MolFromSmiles(leaf["smiles"])
                    if mol:
                        stock_lookup[Chem.MolToSmiles(mol)] = leaf.get("in_stock", False)
                except Exception:
                    pass

            # Render each step in forward order (reactants → product)
            # using expanded_smiles as the reactant list so co-reactants are never dropped
            for si, step in enumerate(reversed(steps_history)):
                target_smi = step.get("target_smiles", "")
                expanded = step.get("expanded_smiles", "")
                cond = _dedup_conditions(step.get("reaction_condition", []))
                rxn_label = classify_reaction(step.get("reaction_template", ""),
                                              step.get("reaction_condition", []),
                                              step.get("reaction_type"))
                rxn_cls = _rxn_type_class(rxn_label)
                is_last_step = (si == len(steps_history) - 1)

                reactants = [r.strip() for r in expanded.split(".") if r.strip()]

                # Step label
                step_num = si + 1
                chain_parts.append(f'<div class="step-row">')
                chain_parts.append(f'<div class="step-label">Step {step_num}</div>')
                chain_parts.append('<div class="step">')

                # Reactant group
                chain_parts.append('<div class="reactants">')
                for ri, reactant in enumerate(reactants):
                    # Check stock status via canonical SMILES lookup
                    stock = None
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            stock = stock_lookup.get(Chem.MolToSmiles(mol))
                    except Exception:
                        pass
                    r_svg = smiles_to_svg_full(reactant, 180, 140)
                    chain_parts.append(_make_mol_box(r_svg, reactant, "", stock=stock))
                    if ri < len(reactants) - 1:
                        chain_parts.append('<span class="plus">+</span>')
                chain_parts.append('</div>')

                # Arrow + conditions
                chain_parts.append(_make_arrow_section(cond, rxn_label, rxn_cls))

                # Product
                tgt_svg = smiles_to_svg_full(target_smi, 200, 160)
                if is_last_step:
                    chain_parts.append(_make_target_mol_box(tgt_svg, target_smi, target_name, str(target_mw)))
                else:
                    chain_parts.append(_make_mol_box(tgt_svg, target_smi))

                chain_parts.append('</div>')  # close .step

                chain_parts.append('</div>')  # close .step-row

            badge_text = "Best Route" if is_best else f"Route {rank}"
            parts.append(f'<div class="route-card {card_cls}">')
            parts.append(f'  <div class="route-header">')
            parts.append(f'    <span class="route-badge">{badge_text}</span>')
            parts.append(f'    <span class="route-score">Score: {score:.4f} · {step_count} step(s)</span>')
            if is_best:
                parts.append(f'    <span class="best-tag">Best</span>')
            parts.append(f'  </div>')
            parts.append(f'  <div class="step-chain">')
            parts.append(f'    {" ".join(chain_parts)}')
            parts.append(f'  </div>')
            parts.append(f'</div>')
            shown += 1

        if shown == 0:
            parts.append(f'<p style="color:#888;text-align:center;margin-top:24px;">No viable route found.</p>')

    parts.append('<div class="footer">This result is a heuristic suggestion, not an experimentally validated protocol. SimpRetro Retrosynthesis Engine</div>')
    parts.append('</body></html>')

    html = "\n".join(parts)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Visualization saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="SimpRetro result visualizer")
    parser.add_argument("input", help="Path to retro_result.json")
    parser.add_argument("-o", "--output", default=None, help="Output HTML path (default: <input>_view.html)")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: file not found: {args.input}")
        return 1

    with open(args.input, "r", encoding="utf-8") as f:
        data = json.load(f)

    output = args.output or os.path.splitext(args.input)[0] + "_view.html"
    generate_flowchart_html(data, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
