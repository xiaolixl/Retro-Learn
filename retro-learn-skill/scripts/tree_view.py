#!/usr/bin/env python3
"""
Tree-style retrosynthesis visualization — top-down branching layout.
Standalone script, does not modify existing visualize.py.

Usage:
  python tree_view.py retro_result.json                    # output: retro_result_tree.html
  python tree_view.py retro_result.json -o custom.html
"""

import argparse
import csv
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# ── Molecule SVG ──────────────────────────────────────────
def smiles_to_svg(smiles, width=160, height=120):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f'<text class="err">?</text>'
    try:
        AllChem.Compute2DCoords(mol)
    except Exception:
        pass
    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    svg = svg.replace('<?xml version="1.0"?>', '')
    svg = svg.replace('<svg', '<svg preserveAspectRatio="xMidYMid meet"')
    return svg

def esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;").replace('"',"&quot;")

# ── Common names ──────────────────────────────────────────
_COMMON_NAMES = {}
def _load_names():
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "common_names.csv")
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8") as f:
            for row in csv.reader(f):
                if len(row) >= 2:
                    _COMMON_NAMES[row[0].strip()] = row[1].strip()
_load_names()

def common_name(smiles):
    if smiles in _COMMON_NAMES:
        return _COMMON_NAMES[smiles]
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            cansmi = Chem.MolToSmiles(mol)
            return _COMMON_NAMES.get(cansmi, "")
    except:
        pass
    return ""

def dedup_conds(clist):
    if not clist:
        return ""
    seen = set()
    u = []
    for c in clist:
        c = str(c).strip()
        if c and c not in seen:
            seen.add(c)
            u.append(c)
    return ", ".join(u)


def classify_reaction(template_str, condition_list, reaction_type=None):
    """Return reaction type label.

    Prefers pre-classified reaction_type from template cache (injected by
    run_retro.py). Falls back to empty string if not classified.
    """
    if reaction_type:
        return reaction_type
    return ""

# ── Layout engine ─────────────────────────────────────────
class TreeLayout:
    """Compute x,y positions for a top-down tree."""

    def __init__(self):
        self.nodes = []       # [(id, x, y, w, h)]
        self.edges = []       # [(from_id, to_id, label)]
        self.counter = 0
        self.MOL_W = 180
        self.MOL_H = 130
        self.GAP_X = 50
        self.GAP_Y = 110
        self.LEAF_GAP = 160
        self.AREA_W = 0
        self.AREA_H = 0

    def new_id(self):
        self.counter += 1
        return f"n{self.counter}"

    def add_node(self, node_id, x, y, w, h):
        self.nodes.append((node_id, x, y, w, h))

    def add_edge(self, frm, to, label=""):
        self.edges.append((frm, to, label))

    def layout_single(self, target, routes):
        """Single-step: target center top, routes below in a row.
        Spacing is computed dynamically from the actual number of routes."""
        tw = self.MOL_W
        th = self.MOL_H
        tx = 0
        ty = 30

        tid = self.new_id()
        self.add_node(tid, tx, ty, tw, th)

        n = len(routes)
        col_gap = max(20, self.GAP_X - (n - 3) * 12) if n > 3 else self.GAP_X
        row_w = n * tw + (n - 1) * col_gap
        start_x = -(row_w / 2) + tw / 2

        for i, route in enumerate(routes):
            rid = self.new_id()
            rx = start_x + i * (tw + col_gap)
            ry = ty + th + self.GAP_Y
            self.add_node(rid, rx, ry, tw, th)

            cond = dedup_conds(route.get("reaction_condition", []))
            self.add_edge(tid, rid, cond)

        self.AREA_W = max(680, row_w + tw * 2)
        self.AREA_H = ty + th + self.GAP_Y + th + 60

        return tid, routes

    def layout_multi(self, all_routes, target):
        """Multi-step: tight columns for intermediates, leaves fan out."""
        tw = self.MOL_W
        th = self.MOL_H
        n_routes = min(len(all_routes), 10)

        # Tight column gap for intermediates
        col_gap = max(25, self.GAP_X - (n_routes - 3) * 15)
        row_w = n_routes * tw + (n_routes - 1) * col_gap
        start_x = -(row_w / 2) + tw / 2

        tid = self.new_id()
        self.add_node(tid, 0, 30, tw, th)

        # Compute leaf spans per route
        leaf_info = []
        for route in all_routes[:n_routes]:
            nl = len(route.get("leaf_reactants", []))
            lg = max(15, col_gap - 20) if nl > 2 else col_gap
            span = nl * tw + max(0, nl - 1) * lg
            leaf_info.append((nl, lg, span))

        # Tight column x positions
        col_x = [start_x + i * (tw + col_gap) for i in range(n_routes)]

        # Iterative push: if leaf blocks overlap, shift right column
        for _pass in range(4):
            moved = False
            for i in range(n_routes - 1):
                _, _, span_i = leaf_info[i]
                _, _, span_j = leaf_info[i + 1]
                right_i = col_x[i] + span_i / 2
                left_j = col_x[i + 1] - span_j / 2
                if right_i + 20 > left_j:
                    col_x[i + 1] += (right_i + 20 - left_j)
                    moved = True
            if not moved:
                break

        # Re-center columns around x=0
        x_shift = -sum(col_x) / n_routes if n_routes else 0
        col_x = [cx + x_shift for cx in col_x]

        route_infos = []
        min_x_all = 0
        max_x_all = 0

        for ri, route in enumerate(all_routes[:n_routes]):
            cx = col_x[ri]
            steps = route.get("steps_history", [])
            leaves = route.get("leaf_reactants", [])
            n_leaves, leaf_gap, _span = leaf_info[ri]

            if not steps:
                continue

            prev_id = tid
            cy = 30 + th + self.GAP_Y

            for si, step in enumerate(steps):
                sid = self.new_id()
                self.add_node(sid, cx, cy, tw, th)
                cond = dedup_conds(step.get("reaction_condition", []))
                self.add_edge(prev_id, sid, cond)
                prev_id = sid
                cy += th + self.GAP_Y
                last_cond = cond  # remember for leaf edge

            # Leaf nodes only needed for multi-leaf routes; single-leaf routes
            # already have the reactant shown in the last step node.
            need_leaves = (n_leaves > 1)
            if not need_leaves:
                pass  # step node itself shows the reactant(s); no extra leaf nodes
            elif n_leaves > 1:
                cy += self.GAP_Y

            leaf_row_w = n_leaves * tw + max(0, n_leaves - 1) * leaf_gap
            lx0 = cx - leaf_row_w / 2 + tw / 2

            for li, leaf in enumerate(leaves):
                if not need_leaves:
                    continue
                lid = self.new_id()
                lx = lx0 + li * (tw + leaf_gap) if n_leaves > 1 else lx0
                self.add_node(lid, lx, cy, tw, th)
                self.add_edge(prev_id, lid, last_cond if last_cond else leaf.get("smiles", "")[:30])
                min_x_all = min(min_x_all, lx)
                max_x_all = max(max_x_all, lx + tw)

            route_infos.append({
                "rank": route.get("route_rank", ri + 1),
                "score": route.get("route_score", 0),
                "end_y": cy + th,
            })

        max_y = max((r["end_y"] for r in route_infos), default=0) + 60
        self.AREA_W = max(680, (max_x_all - min_x_all) + tw * 2 + 80)
        self.AREA_H = max(max_y + 80, 600)
        return tid, route_infos

# ── SVG render ────────────────────────────────────────────
def render_svg(layout, target_smiles, routes_data, target_svg, mode):
    """Generate the full SVG tree."""
    W = max(layout.AREA_W, 680)
    H = layout.AREA_H
    cx = W / 2

    # Collect node metadata: (node_id, smiles, name, stock, is_target)
    node_meta = {}

    def t(x):
        return x + cx  # translate to center

    parts = [f'''<svg viewBox="0 0 {W} {H}" width="{W}" role="img"
  xmlns="http://www.w3.org/2000/svg">''']

    # Build node metadata
    for nid, x, y, w, h in layout.nodes:
        node_meta[nid] = {"x": t(x), "y": y, "w": w, "h": h, "smiles": "", "name": "", "stock": None}

    # Fill in SMILES from routes
    # First node is target
    target_nid = layout.nodes[0][0]
    node_meta[target_nid]["smiles"] = target_smiles
    node_meta[target_nid]["is_target"] = True

    # Map route nodes
    node_idx = 1
    if mode == "single_step":
        for ri, route in enumerate(routes_data):
            if node_idx >= len(layout.nodes):
                break
            nid = layout.nodes[node_idx][0]
            reactants = route.get("reactants", [])
            if len(reactants) == 1:
                node_meta[nid]["smiles"] = reactants[0]["smiles"]
                node_meta[nid]["stock"] = reactants[0].get("in_stock")
            elif len(reactants) > 1:
                node_meta[nid]["smiles"] = " + ".join(r["smiles"] for r in reactants)
                node_meta[nid]["multipart"] = [(r["smiles"], r.get("in_stock")) for r in reactants]
            node_meta[nid]["score"] = route.get("score", 0)
            node_meta[nid]["rxn_type"] = classify_reaction(
                route.get("reaction_template", ""),
                route.get("reaction_condition", []),
                route.get("reaction_type"))
            node_idx += 1

    else:
        best_score = max((r.get("route_score", -999) for r in routes_data), default=0)
        best_count = sum(1 for r in routes_data if r.get("route_score", -999) == best_score)
        only_one_best = (best_count == 1)
        route_idx = 0
        for route in routes_data[:10]:
            route_idx += 1
            is_best = only_one_best and route.get("route_score", -999) == best_score
            steps = route.get("steps_history", [])
            leaves = route.get("leaf_reactants", [])
            n_leaves = len(leaves)
            n_steps = len(steps)
            need_leaves = (n_leaves > 1)

            for si, step in enumerate(steps):
                if node_idx >= len(layout.nodes):
                    break
                nid = layout.nodes[node_idx][0]
                expanded = step.get("expanded_smiles", "")
                is_last = (si == n_steps - 1)
                if "." in expanded:
                    if is_last and not need_leaves:
                        # Last step, no leaf nodes: render as multipart directly
                        parts_smi = expanded.split(".")
                        node_meta[nid]["multipart"] = [(s, False) for s in parts_smi]
                    else:
                        # Dot-separated but not last/no-multipart: compact text
                        node_meta[nid]["smiles"] = expanded.replace(".", " · ")
                else:
                    node_meta[nid]["smiles"] = expanded
                node_meta[nid]["score"] = step.get("step_score", 0)
                node_meta[nid]["route_score"] = route.get("route_score", 0)
                node_meta[nid]["is_best"] = is_best
                node_meta[nid]["route_label"] = f"Route {route_idx}"
                node_meta[nid]["rxn_type"] = classify_reaction(
                    step.get("reaction_template", ""),
                    step.get("reaction_condition", []),
                    step.get("reaction_type"))
                node_idx += 1

            for leaf in leaves:
                if not need_leaves:
                    continue
                if node_idx >= len(layout.nodes):
                    break
                nid = layout.nodes[node_idx][0]
                node_meta[nid]["smiles"] = leaf["smiles"]
                node_meta[nid]["stock"] = leaf.get("in_stock")
                node_meta[nid]["route_score"] = route.get("route_score", 0)
                node_meta[nid]["is_best"] = is_best
                node_meta[nid]["route_label"] = f"Route {route_idx}"
                node_idx += 1

    # Draw edges
    for frm, to, label in layout.edges:
        if frm not in node_meta or to not in node_meta:
            continue
        fm = node_meta[frm]
        tm = node_meta[to]
        x1 = fm["x"] + fm["w"] / 2
        y1 = fm["y"] + fm["h"]
        x2 = tm["x"] + tm["w"] / 2
        y2 = tm["y"]

        # Draw vertical + horizontal connector
        mid_y = (y1 + y2) / 2
        color = "#534AB7"

        parts.append(f'<polyline points="{x1},{y1} {x1},{mid_y} {x2},{mid_y} {x2},{y2}"')
        parts.append(f'  fill="none" stroke="{color}" stroke-width="1.5"/>')

        # Condition label: right of vertical drop, except leftmost route goes left
        if label:
            ly = (mid_y + y2) / 2 - 4
            on_left = x2 < cx  # left of tree center
            if on_left:
                parts.append(f'<text x="{x2 - 6}" y="{ly}" class="edge-label" text-anchor="end">{esc(label)[:55]}</text>')
            else:
                parts.append(f'<text x="{x2 + 6}" y="{ly}" class="edge-label">{esc(label)[:55]}</text>')

    # Draw nodes
    for nid, meta in node_meta.items():
        x, y, w, h = meta["x"], meta["y"], meta["w"], meta["h"]
        is_target = meta.get("is_target", False)
        smiles = meta.get("smiles", "")
        stock = meta.get("stock")
        score = meta.get("score")
        multipart = meta.get("multipart")

        # Card background
        fill = "#EEEDFE" if is_target else "#FFFFFF"
        stroke = "#534AB7" if is_target else "#B4B2A9"
        if stock is True:
            fill = "#E1F5EE"
            stroke = "#0F6E56"

        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="10"')
        parts.append(f'  fill="{fill}" stroke="{stroke}" stroke-width="1"/>')

        if multipart:
            # Multiple molecules in one box
            subw = (w - 16) // len(multipart)
            for mi, (sm, st) in enumerate(multipart):
                sx = x + 8 + mi * subw
                mol_svg = smiles_to_svg(sm, subw, h - 40)
                mol_svg = mol_svg.replace('<svg', f'<svg x="{sx}" y="{y + 2}" width="{subw}" height="{h - 40}"')
                parts.append(mol_svg)
        else:
            # Single molecule
            mol_svg = smiles_to_svg(smiles, w - 8, h - 34)
            mol_svg = mol_svg.replace('<svg', f'<svg x="{x + 4}" y="{y + 2}" width="{w - 8}" height="{h - 34}"')
            parts.append(mol_svg)

        # Label below structure
        ly = y + h - 22
        cname = common_name(smiles)
        display = cname if cname else smiles
        parts.append(f'<text x="{x + w / 2}" y="{ly}" class="node-label" text-anchor="middle">{esc(display)[:28]}</text>')

        # Reaction type tag
        rxn = meta.get("rxn_type", "")
        if rxn:
            parts.append(f'<text x="{x + w / 2}" y="{y + h - 6}" class="rxn-label" text-anchor="middle">{esc(rxn)[:30]}</text>')

        # Route label (above first node of each route)
        route_label = meta.get("route_label", "")
        if route_label and not is_target:
            parts.append(f'<rect x="{x + 2}" y="{y - 14}" width="52" height="14" rx="3" fill="#534AB7"/>')
            parts.append(f'<text x="{x + 28}" y="{y - 4}" class="stock-badge" text-anchor="middle">{esc(route_label)}</text>')

        # Target badge (above node)
        if is_target:
            parts.append(f'<text x="{x + w / 2}" y="{y - 8}" class="target-label" text-anchor="middle">Target Molecule</text>')

        # Stock badge
        if stock is True:
            parts.append(f'<rect x="{x + w - 40}" y="{y + 2}" width="36" height="14" rx="3" fill="#0F6E56"/>')
            parts.append(f'<text x="{x + w - 22}" y="{y + 9}" class="stock-badge" text-anchor="middle">In Stock</text>')

        if score is not None and not is_target:
            parts.append(f'<text x="{x + 4}" y="{y + 12}" class="score-label">{score:.3f}</text>')

        # Best route badge (multi-step)
        if meta.get("is_best"):
            parts.append(f'<rect x="{x + w - 48}" y="{y + 2}" width="44" height="14" rx="3" fill="#0F6E56"/>')
            parts.append(f'<text x="{x + w - 26}" y="{y + 9}" class="stock-badge" text-anchor="middle">Best</text>')

    parts.append('</svg>')
    return "\n".join(parts)

# ── HTML wrapper ──────────────────────────────────────────
HTML_TPL = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Retrosynthesis Tree View</title>
<style>
  body { font-family: -apple-system, 'Segoe UI', sans-serif; background: #f8f8f6; color: #2c2c2a; padding: 24px; margin: 0; text-align: center; }
  h1 { font-size: 18px; font-weight: 500; color: #26215C; margin-bottom: 8px; }
  .subtitle { font-size: 13px; color: #888; margin-bottom: 20px; }
  .svg-wrap { max-width: 100%; margin: 0 auto; overflow-x: auto; overflow-y: hidden; padding: 10px; -webkit-overflow-scrolling: touch; }
  svg { font-family: -apple-system, 'Segoe UI', sans-serif; }
  .node-label { font-size: 11px; fill: #2c2c2a; }
  .rxn-label { font-size: 9px; fill: #534AB7; font-weight: 500; }
  .target-label { font-size: 10px; fill: #534AB7; font-weight: 500; }
  .edge-label { font-size: 9px; fill: #534AB7; }
  .stock-badge { font-size: 8px; fill: #fff; }
  .score-label { font-size: 9px; fill: #888; }
  .err { font-size: 10px; fill: #999; }
  .footer { font-size: 11px; color: #999; margin-top: 24px; border-top: 1px solid #eee; padding-top: 12px; }
</style>
</head>
<body>
<h1>Retrosynthesis Tree View</h1>
<div class="subtitle">{subtitle}</div>
<div class="svg-wrap">
{svg_content}
</div>
<div class="footer">SimpRetro Retrosynthesis Engine — heuristic suggestions, not experimentally validated.</div>
</body>
</html>"""

# ── Main ──────────────────────────────────────────────────
def generate(data, output_path):
    target = data["data"]["target_molecule"]
    target_smiles = target["smiles"]
    mode = data["data"].get("mode", "single_step")

    layout = TreeLayout()

    if mode == "single_step":
        routes = data["data"].get("retrosynthesis_routes", [])
        layout.layout_single(target, routes)
        subtitle = f"Target: {target_smiles} · {len(routes)} single-step route(s)"
        svg = render_svg(layout, target_smiles, routes, None, "single_step")

    else:
        all_routes = data["data"].get("all_routes", [])
        recommended = data["data"].get("recommended_route")
        if not all_routes and recommended:
            all_routes = [{"route_rank": 1, "route_score": recommended.get("route_score", 0),
                           "steps": recommended.get("actual_steps", 0),
                           "leaf_reactants": recommended.get("leaf_reactants", []),
                           "steps_history": recommended.get("steps", [])}]
        layout.layout_multi(all_routes, target)
        n = len(all_routes)
        subtitle = f"Target: {target_smiles} · {n} multi-step route(s)"
        svg = render_svg(layout, target_smiles, all_routes, None, "multi_step")

    html = HTML_TPL.replace("{subtitle}", subtitle).replace("{svg_content}", svg)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Tree view saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Retrosynthesis tree visualizer")
    parser.add_argument("input", help="Path to retro_result.json")
    parser.add_argument("-o", "--output", default=None, help="Output HTML path")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: file not found: {args.input}")
        return 1

    with open(args.input, "r", encoding="utf-8") as f:
        data = json.load(f)

    output = args.output or os.path.splitext(args.input)[0] + "_tree.html"
    generate(data, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
