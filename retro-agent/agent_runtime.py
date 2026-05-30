import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import uuid4

from openai import OpenAI
from pydantic import BaseModel, Field

from chem_resolution import resolve_identifier, resolve_identifier_list

# Import engine from retro-agent/engine/ (local copy, independent of retro-learn-skill)
_AGENT_ROOT = Path(__file__).resolve().parent
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))
from engine.retro_engine import DEFAULT_DATABASE, DEFAULT_WEIGHTS, run_retrosynthesis
from engine.route_planner import plan_retrosynthesis


PROJECT_ROOT = Path(__file__).resolve().parent
ENGINE_ROOT = PROJECT_ROOT.parent / "SimpRetro4Learn"  # submodule — data files only
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "user_output" / "agent_runs"
CACHE_ROOT = PROJECT_ROOT / "user_output" / "cached_results"
DEFAULT_MODEL = os.getenv("OPENAI_MODEL", "gpt-5.4-mini")


class ParsedAgentRequest(BaseModel):
    language: str = "zh"
    step_count: int = 1
    top_k: Optional[int] = None
    target_smiles: Optional[str] = None
    target_name: Optional[str] = None
    preferred_reactants: List[str] = Field(default_factory=list)
    preferred_reactant_names: List[str] = Field(default_factory=list)
    database: Optional[str] = None
    weights: Optional[List[float]] = None
    needs_clarification: bool = False
    clarification_question: Optional[str] = None
    summary_of_request: str = ""


class ResolvedAgentRequest(BaseModel):
    language: str = "zh"
    step_count: int = 1
    top_k: int = 5
    target_smiles: str
    preferred_reactants: List[str] = Field(default_factory=list)
    database: str = DEFAULT_DATABASE
    weights: List[float] = Field(default_factory=lambda: list(DEFAULT_WEIGHTS))
    summary_of_request: str = ""
    resolution_notes: List[str] = Field(default_factory=list)


def _extract_json_object(text: str) -> Dict[str, Any]:
    stripped = text.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```(?:json)?", "", stripped).strip()
        stripped = re.sub(r"```$", "", stripped).strip()
    try:
        return json.loads(stripped)
    except json.JSONDecodeError:
        match = re.search(r"\{.*\}", stripped, flags=re.DOTALL)
        if not match:
            raise ValueError(f"Could not parse JSON from model output: {text}")
        return json.loads(match.group(0))


def _compact_planning_payload(result: Dict[str, Any]) -> Dict[str, Any]:
    data = result["data"]
    payload = {
        "mode": data["mode"],
        "target_molecule": data["target_molecule"],
        "query_constraints": data.get("query_constraints", data.get("request_params", {})),
    }
    if data["mode"] == "single_step":
        routes = data.get("recommended_routes") or data.get("retrosynthesis_routes", [])
        payload["recommended_routes"] = [
            {
                "route_rank": route.get("route_rank", i + 1),
                "score": route.get("score", 0),
                "all_reactants_in_stock": route.get("all_reactants_in_stock", False),
                "matched_preferred_reactants": route.get("matched_preferred_reactants", []),
                "reaction_condition": route.get("reaction_condition", []),
                "reactants": route.get("reactants", []),
            }
            for i, route in enumerate(routes)
        ]
    else:
        route = data.get("recommended_route")
        if route is None:
            payload["recommended_route"] = None
        else:
            payload["recommended_route"] = {
                "requested_steps": route.get("requested_steps", route.get("actual_steps", 0)),
                "actual_steps": route.get("actual_steps", len(route.get("steps", []))),
                "completed_requested_steps": route.get("completed_requested_steps", True),
                "route_score": route.get("route_score", 0),
                "average_step_score": route.get("average_step_score", 0),
                "stock_leaf_ratio": route.get("stock_leaf_ratio", 0),
                "matched_preferred_reactants": route.get("matched_preferred_reactants", []),
                "final_leaf_reactants": route.get("leaf_reactants") or route.get("final_leaf_reactants", []),
                "steps": [
                    {
                        "step_number": step.get("step_number", i + 1),
                        "expanded_smiles": step.get("expanded_smiles", ""),
                        "step_score": step.get("step_score", 0),
                        "matched_preferred_reactants": step.get("matched_preferred_reactants", []),
                        "reaction_condition": step.get("reaction_condition", []),
                        "reactants": [{"smiles": s} for s in step.get("expanded_smiles", "").split(".")],
                    }
                    for i, step in enumerate(route.get("steps", []))
                ],
            }
    return payload


class RetrosynthesisAgent:
    def __init__(
        self,
        model: Optional[str] = None,
        output_root: Path = DEFAULT_OUTPUT_ROOT,
    ) -> None:
        api_key = os.getenv("OPENAI_API_KEY")
        if not api_key:
            raise ValueError("OPENAI_API_KEY is not set.")

        base_url = os.getenv("OPENAI_BASE_URL")
        if base_url:
            self.client = OpenAI(api_key=api_key, base_url=base_url)
        else:
            self.client = OpenAI(api_key=api_key)

        self.model = model or DEFAULT_MODEL
        self.output_root = output_root
        self.base_url = base_url or ""

        # Load template→reaction_type cache (populated by LLM at runtime, grows over time)
        self.template_type_map: Dict[str, str] = {}
        self._template_cache_path = PROJECT_ROOT / "template_reaction_types.json"
        self._template_cache_dirty = False
        if self._template_cache_path.exists():
            try:
                self.template_type_map = json.loads(self._template_cache_path.read_text(encoding="utf-8"))
            except Exception:
                pass

    def _uses_deepseek_chat_api(self) -> bool:
        model_name = (self.model or "").lower()
        base_url = self.base_url.lower()
        return "deepseek" in model_name or "deepseek.com" in base_url

    def _llm_text(self, prompt: str) -> str:
        if self._uses_deepseek_chat_api():
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "user", "content": prompt},
                ],
                stream=False,
            )
            return response.choices[0].message.content or ""

        response = self.client.responses.create(
            model=self.model,
            input=prompt,
        )
        return response.output_text

    def _llm_json(self, prompt: str) -> Dict[str, Any]:
        if self._uses_deepseek_chat_api():
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": "Return valid json only.",
                    },
                    {
                        "role": "user",
                        "content": prompt,
                    },
                ],
                response_format={"type": "json_object"},
                stream=False,
            )
            content = response.choices[0].message.content or "{}"
            return _extract_json_object(content)

        return _extract_json_object(self._llm_text(prompt))

    def parse_user_request(self, user_message: str) -> ParsedAgentRequest:
        prompt = f"""
You are a chemistry agent request parser.
Extract a single JSON object from the user's request for a retrosynthesis assistant.

Scenario detection rules (check the user's wording):
- **Default / Scenario A**: If the user only provides a molecule name or SMILES (with no explicit step count, starting material, or carbon constraint), treat it as a SINGLE-STEP retrosynthesis. Set step_count=1 and preferred_reactants=[]. Examples: "retrosynthesis for acetone", "CC(=O)C", "analyze stilbene". This is the most common case.
- **Scenario B**: If the user explicitly asks for multi-step (e.g. "3-step", "multi-step", "多步") with NO starting material and NO carbon constraint, set step_count=N and preferred_reactants=[].
- **Scenario C**: If the user specifies a starting material (e.g. "from toluene", "starting from benzene", "合成…从…") OR a carbon-count constraint (e.g. "≤4C", "不超过4碳", "from 3-carbon building blocks"), set step_count=5 (default unless user specifies) and put the starting material SMILES in preferred_reactants. This is a multi-step retrosynthesis with constraints.
  - **When starting material is specified WITHOUT carbon constraint**: set database to "emol_under_0_carbons" (unrestricted). The starting material acts as the stock constraint — the engine will prioritize routes that reach it.
  - **When carbon constraint is specified WITHOUT starting material**: set database to the matching emol DB below.
  - **When BOTH starting material AND carbon constraint are specified**: set database to the matching emol DB. Both constraints apply.

- **Scenario B** (multi-step, NO starting material, NO carbon constraint): set step_count=3, database="emol_under_0_carbons" (unrestricted). Multi-step search needs the full stock DB to find viable routes.

Carbon-count constraint mapping (only applies when user explicitly specifies a carbon limit):
- Detect phrases like "≤4C", "4碳以内", "不超过4碳", "from ≤4 carbon", "under 4 carbons", "三碳以内" etc.
- Map to the `database` field using this table:

| Max carbons | database value              |
|-------------|---------------------------|
| ≤3C         | emol_under_3_carbons      |
| ≤4C         | emol_under_4_carbons      |
| ≤5C         | emol_under_5_carbons      |
| ≤6C         | emol_under_6_carbons      |
| Unrestricted / 不限制 | emol_under_0_carbons |

- When a carbon constraint is detected, also set step_count=3 (multi-step) unless the user explicitly says single-step.

Top-K extraction:
- If the user specifies "top N", "top-N", "top_k=N" or similar (e.g. "top 5", "top-3"), set top_k=N.
- If the user does not specify, set top_k=null (engine default applies).

Name-to-SMILES resolution:
- If the user provides a molecule name (e.g. "stilbene", "toluene", "aspirin", "acetone"), convert it to SMILES using your chemical knowledge. Put the SMILES in target_smiles (for target) or preferred_reactants (for starting materials). Also record the original name in target_name / preferred_reactant_names.
- **Stereochemistry preservation (CRITICAL)**: Use ISOMERIC SMILES whenever the molecule has known stereochemistry. Use your organic chemistry knowledge to decide:
  - If the molecule name describes a specific isomer (e.g. "trans-stilbene", "L-alanine", "(R)-ibuprofen") → include / \ @ @@ markers.
  - If the common name conventionally refers to a specific stereoisomer in chemistry (e.g. "stilbene" is trans/E by default, "glucose" is D-glucose) → include the markers.
  - If the molecule genuinely has NO stereochemistry (e.g. acetone, toluene, aspirin, benzene) → use non-isomeric SMILES.
  - Do NOT fabricate stereochemistry for achiral molecules.
- If the user gives an explicit SMILES string, use it directly.
- Only set needs_clarification=true if you truly cannot identify the molecule.

Other:
- Return JSON only.
- If the request language is Chinese, set language to "zh"; otherwise use "en".
- weights should be null unless the user explicitly specifies them.

JSON schema:
{{
  "language": "zh or en",
  "step_count": 1,
  "top_k": "integer or null",
  "target_smiles": "string (canonical SMILES) or null",
  "target_name": "original name if user gave a name, or null",
  "preferred_reactants": ["canonical SMILES", "..."],
  "preferred_reactant_names": ["original name", "..."],
  "database": "string or null (e.g. emol_under_4_carbons)",
  "weights": [0.1, 0.2, 0.5, 0.0] or null,
  "needs_clarification": false,
  "clarification_question": "string or null",
  "summary_of_request": "short summary"
}}

User request:
{user_message}
""".strip()
        parsed = self._llm_json(prompt)
        return ParsedAgentRequest.model_validate(parsed)

    def resolve_request(self, parsed_request: ParsedAgentRequest) -> ResolvedAgentRequest:
        if parsed_request.needs_clarification:
            question = parsed_request.clarification_question or "Please provide the target molecule as a SMILES string."
            raise ValueError(question)

        # LLM now resolves names to SMILES directly; use its output first
        if parsed_request.target_smiles:
            target_smiles = parsed_request.target_smiles
            target_source = "llm"
        elif parsed_request.target_name:
            target_smiles, target_source = resolve_identifier(parsed_request.target_name)
            if target_smiles is None:
                raise ValueError(
                    "I could not reliably resolve the target molecule to a structure. Please provide the target molecule as a SMILES string."
                )
        else:
            raise ValueError("Please provide the target molecule as a SMILES string or a resolvable molecule name.")

        # ── RDKit canonicalization with stereochemistry preservation ──
        # Canonicalize to catch any SMILES issues; use isomericSmiles=True to
        # preserve E/Z and chiral markers that the LLM included.
        from rdkit import Chem
        mol = Chem.MolFromSmiles(target_smiles)
        if mol is not None:
            target_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

        preferred_inputs = list(parsed_request.preferred_reactants)
        preferred_reactants, unresolved_preferred = resolve_identifier_list(preferred_inputs)
        if unresolved_preferred:
            unresolved_display = ", ".join(unresolved_preferred)
            raise ValueError(
                f"I could not reliably resolve these preferred reactants: {unresolved_display}. Please provide their SMILES strings."
            )

        resolution_notes = [f"target_source={target_source}"]
        return ResolvedAgentRequest(
            language=parsed_request.language or "zh",
            step_count=max(1, parsed_request.step_count),
            top_k=parsed_request.top_k or 5,
            target_smiles=target_smiles,
            preferred_reactants=preferred_reactants,
            database=parsed_request.database or DEFAULT_DATABASE,
            weights=parsed_request.weights or list(DEFAULT_WEIGHTS),
            summary_of_request=parsed_request.summary_of_request,
            resolution_notes=resolution_notes,
        )

    def _prepare_output_dir(self) -> Path:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_id = uuid4().hex[:8]
        output_dir = self.output_root / f"{timestamp}_{run_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def explain_results(self, user_message: str, resolved_request: ResolvedAgentRequest, planning_result: Dict[str, Any]) -> str:
        payload = _compact_planning_payload(planning_result)
        lang_dir = "Chinese" if resolved_request.language == "zh" else "English"
        # Only pass essential request info, not full dump
        req_summary = {
            "target_smiles": resolved_request.target_smiles,
            "step_count": resolved_request.step_count,
            "preferred_reactants": resolved_request.preferred_reactants,
            "database": resolved_request.database,
            "summary": resolved_request.summary_of_request,
        }
        prompt = f"""
You are a helpful retrosynthesis agent.
Write the answer in {lang_dir}.

## Output Structure

### 1. Target Information
Start with the target molecule:
- Name (if known from common_names or the user's request)
- SMILES string
- Molecular weight (if available in planning_result)

### 2. Result Summary
- If mode is single_step: summarize the top routes in rank order (up to 3). For each route, briefly state the reaction type, reactants, and score.
- If mode is multi_step: present the best route step-by-step. Show each step's reactants → conditions → product, with scores.

### 3. Key Observations
- Preferred reactant matches: if the user specified starting materials, report whether any routes reached them.
- Stock status: mention whether reactants are in-stock or need synthesis.
- If no viable routes were found, say so clearly.

### 4. Disclaimer
End with this exact warning (translated to {lang_dir} if needed):
"⚠️ 以上路线为计算辅助建议，未经实验验证，实际合成前需经专业人员审核。"
Or in English: "⚠️ The routes above are computational suggestions, not experimentally validated. Expert review is required before experimental use."

## Style Rules
- Use markdown headings (##) and bullets for readability.
- Keep the answer practical and concise — avoid repeating raw data verbatim.
- Do NOT fabricate information not present in the planning result.
- Reaction type labels: only name known types (oxidation, reduction, elimination, Diels-Alder, Wittig, Grignard, etc.). If uncertain, omit the type label.

Original user request:
{user_message}

Request summary:
{json.dumps(req_summary, indent=2, ensure_ascii=False)}

Planning result summary:
{json.dumps(payload, indent=2, ensure_ascii=False)}
""".strip()
        return self._llm_text(prompt).strip()

    def _assign_reaction_types(self, planning_result: Dict[str, Any]) -> None:
        """Two-phase reaction type assignment: cache lookup + LLM fallback.

        Phase 1: O(1) template cache lookup (template_reaction_types.json).
        Phase 2: For cache misses, sends concrete SMILES + conditions to LLM
                 in a single batch call. Results are saved back to the cache
                 so future runs skip LLM.

        Mutates planning_result in-place.
        """
        data = planning_result.get("data", {})
        mode = data.get("mode", "single_step")
        target_smiles = data.get("target_molecule", {}).get("smiles", "")

        # ── Collect unique step groups (reactants, product, conditions) ──
        # Each group maps to one or more step dicts that share the same transformation.
        StepGroup = Dict[str, Any]  # {reactants, product, conditions, templates, step_refs}
        groups: Dict[tuple, StepGroup] = {}

        def _add_step(reactants, product, conditions, template, step_ref):
            key = (
                tuple(sorted(reactants)),
                product,
                tuple(sorted(str(c) for c in (conditions or []))),
            )
            if key not in groups:
                groups[key] = {
                    "reactants": reactants,
                    "product": product,
                    "conditions": conditions or [],
                    "templates": set(),
                    "step_refs": [],
                }
            groups[key]["templates"].add(template.strip() if template else "")
            groups[key]["step_refs"].append(step_ref)

        if mode == "single_step":
            for route in data.get("retrosynthesis_routes", []):
                _add_step(
                    [r["smiles"] for r in route.get("reactants", [])],
                    target_smiles,
                    route.get("reaction_condition", []),
                    route.get("reaction_template", ""),
                    route,
                )
        else:
            rec = data.get("recommended_route")
            if rec:
                for step in rec.get("steps", []):
                    _add_step(
                        [s.strip() for s in step.get("expanded_smiles", "").split(".") if s.strip()],
                        step.get("target_smiles", ""),
                        step.get("reaction_condition", []),
                        step.get("reaction_template", ""),
                        step,
                    )
            for route in data.get("all_routes", []):
                for step in route.get("steps_history", []):
                    _add_step(
                        [s.strip() for s in step.get("expanded_smiles", "").split(".") if s.strip()],
                        step.get("target_smiles", ""),
                        step.get("reaction_condition", []),
                        step.get("reaction_template", ""),
                        step,
                    )

        if not groups:
            return

        # ── Phase 1: Cache lookup ──
        missed_groups = []  # list of (group_key, group)
        for key, group in groups.items():
            # Try each template in this group against the cache
            cached_type = ""
            for tmpl in group["templates"]:
                if tmpl and tmpl in self.template_type_map:
                    ct = self.template_type_map[tmpl]
                    if ct and ct != "Other":
                        cached_type = ct
                        break
            if cached_type:
                for ref in group["step_refs"]:
                    ref["reaction_type"] = cached_type
            else:
                missed_groups.append((key, group))

        if not missed_groups:
            return  # All cache hits — done

        # ── Phase 2: LLM classify cache misses ──
        group_keys = [k for k, _ in missed_groups]
        steps_text_parts = []
        for i, (key, group) in enumerate(missed_groups):
            r_str = " + ".join(group["reactants"])
            c_str = ", ".join(group["conditions"]) if group["conditions"] else "(none)"
            steps_text_parts.append(f"Step {i}: {r_str}  ->  {group['product']}  |  Conditions: {c_str}")

        prompt = f"""You are an expert organic chemist. Identify the reaction type for each step below.

For each step you are given: Reactants -> Product | Conditions

Classify each reaction concisely. Use specific subtypes when the mechanism is clear:
- "Cross-coupling (Gilman cuprate)" for organocopper + organic halide -> coupled product
- "Diels-Alder cycloaddition" for conjugated diene + dienophile -> cyclohexene derivative
- "Fischer esterification" for acid + alcohol -> ester (+ H2O), esp. with H+
- "Hydrolysis" for ester/amide/nitrile + H2O -> acid + alcohol/amine
- "Oxidation" for gain of O or loss of H
- "Reduction" for gain of H or loss of O
- "Elimination" for C=C formation with loss of H2O/HX
- "Substitution (SN2)" for bimolecular nucleophilic substitution
- "Grignard addition" for RMgX + carbonyl -> alcohol
- "Wittig reaction" for Ph3P=CR2 + carbonyl -> alkene
- "Aldol condensation" for enolate + carbonyl -> unsaturated carbonyl
- "Catalytic hydrogenation" for H2 + C=C with Pd/Pt/Ni
- "Halogenation" for Br2/Cl2 addition across C=C
- "Alkylation" for introduction of alkyl group
- Other specific types as appropriate

Return JSON only:
{{"classifications": [{{"step_id": 0, "reaction_type": "Cross-coupling (Gilman cuprate)"}}, ...]}}

Reaction steps:
{chr(10).join(steps_text_parts)}"""

        try:
            result = self._llm_json(prompt)
            for c in result.get("classifications", []):
                idx = c.get("step_id", -1)
                rxn_type = c.get("reaction_type", "")
                if isinstance(idx, int) and 0 <= idx < len(group_keys):
                    key = group_keys[idx]
                    group = groups[key]
                    # Inject into all steps sharing this group
                    for ref in group["step_refs"]:
                        ref["reaction_type"] = rxn_type
                    # Save to cache for each template in this group
                    for tmpl in group["templates"]:
                        if tmpl and rxn_type:
                            self.template_type_map[tmpl] = rxn_type
                            self._template_cache_dirty = True
        except Exception:
            pass  # LLM call failed — steps remain unclassified, fall back to hard-coded

        # Persist cache if it was updated
        if self._template_cache_dirty:
            try:
                self._template_cache_path.write_text(
                    json.dumps(self.template_type_map, indent=2, ensure_ascii=False),
                    encoding="utf-8",
                )
                self._template_cache_dirty = False
            except Exception:
                pass

    @staticmethod
    def _scenario_c_needs_llm_fallback(planning_result: Dict[str, Any]) -> bool:
        """Check if Scenario C has no route that reached preferred reactants."""
        data = planning_result.get("data", {})
        if data.get("mode") != "multi_step":
            return False
        all_routes = data.get("all_routes", [])
        if not all_routes:
            return True  # No routes at all
        # Check if any route completed (reached preferred reactants)
        # The engine marks completed routes via 'completed' field
        # If no completed flag, check if recommended_route exists
        rec = data.get("recommended_route")
        return rec is None

    def _llm_design_routes(self, target_smiles: str, target_name: str,
                           preferred_smiles_list: list, preferred_names: list,
                           step_count: int, language: str) -> list:
        """LLM-designed synthetic routes when engine fails to reach preferred reactants.

        Returns a list of route dicts compatible with planning_result['data']['all_routes'].
        """
        target_label = target_name or target_smiles
        starting_label = ", ".join(preferred_names) if preferred_names else ", ".join(preferred_smiles_list)
        lang_hint = "Chinese" if language == "zh" else "English"

        prompt = f"""You are an expert synthetic organic chemist. Design 1-2 practical synthetic routes.

Target molecule: {target_label} ({target_smiles})
Starting material(s): {starting_label} ({', '.join(preferred_smiles_list)})
Maximum steps: {step_count}

For each route, describe the forward synthesis step-by-step. For each step, provide:
- Reactants (canonical SMILES)
- Product (canonical SMILES)
- Reaction conditions (reagents, solvent, temperature)
- Reaction type/name

Return JSON only:
{{
  "routes": [
    {{
      "summary": "brief description of the route strategy",
      "steps": [
        {{
          "step_number": 1,
          "reactants": ["Cc1ccccc1"],
          "product": "O=Cc1ccccc1",
          "conditions": "KMnO4, H2SO4, heat",
          "reaction_type": "Oxidation"
        }},
        ...
      ]
    }}
  ]
}}

Write the response in {lang_hint}. Use standard organic chemistry knowledge — avoid exotic or untested transformations.
""".strip()

        try:
            result = self._llm_json(prompt)
            routes = result.get("routes", [])
        except Exception:
            return []

        # Convert LLM output to planning_result route format
        converted = []
        for rank, route in enumerate(routes[:2], start=1):
            steps_data = route.get("steps", [])
            if not steps_data:
                continue

            leaf_reactants = []
            seen_leaves = set()
            forward_steps = []
            all_intermediates = set()

            for i, step in enumerate(steps_data):
                step_num = i + 1
                reactants = step.get("reactants", [])
                product = step.get("product", "")
                conditions_raw = step.get("conditions", "")
                rxn_type = step.get("reaction_type", "")

                # Split conditions string into list
                if isinstance(conditions_raw, str):
                    conditions = [c.strip() for c in conditions_raw.split(",") if c.strip()]
                elif isinstance(conditions_raw, list):
                    conditions = conditions_raw
                else:
                    conditions = []

                # Track which reactants are external (not produced by earlier steps)
                for r in reactants:
                    if r not in all_intermediates and r not in seen_leaves:
                        seen_leaves.add(r)
                        leaf_reactants.append({"smiles": r, "in_stock": True, "molecular_weight": 0})

                all_intermediates.add(product)

                forward_steps.append({
                    "target_smiles": product,
                    "expanded_smiles": ".".join(reactants),
                    "reaction_template": "",
                    "reaction_condition": conditions,
                    "step_score": 0,
                    "reaction_type": rxn_type,
                })

            # Engine stores steps in retrosynthetic order (step 1 = last forward step).
            # visualize.py reverses for display. So we reverse forward_steps and renumber.
            steps_history = list(reversed(forward_steps))
            for i, s in enumerate(steps_history):
                s["step_number"] = i + 1

            converted.append({
                "route_rank": rank + 900,
                "route_score": 0,
                "steps": len(steps_data),
                "source": "LLM-designed",
                "leaf_reactants": leaf_reactants,
                "steps_history": steps_history,
            })

        return converted

    def run(self, user_message: str) -> Dict[str, Any]:
        parsed_request = self.parse_user_request(user_message)
        output_dir = self._prepare_output_dir()
        (output_dir / "parsed_request.json").write_text(
            json.dumps(parsed_request.model_dump(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        try:
            resolved_request = self.resolve_request(parsed_request)
        except ValueError as exc:
            clarification = {
                "status": "needs_clarification",
                "message": str(exc),
                "reply_markdown": str(exc),
                "parsed_request": parsed_request.model_dump(),
                "output_dir": str(output_dir),
            }
            (output_dir / "agent_result.json").write_text(
                json.dumps(clarification, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            return clarification

        (output_dir / "resolved_request.json").write_text(
            json.dumps(resolved_request.model_dump(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        # ── Check cache by target SMILES + scenario key ──
        cache_key = self._cache_key(resolved_request)
        cached = self._load_cache(cache_key)
        if cached is not None:
            # Copy cached files into this run's output dir for frontend serving
            for fname in ("planning_result.json", "route_view.html", "agent_reply.md"):
                src = CACHE_ROOT / cache_key / fname
                if src.exists():
                    shutil.copy2(str(src), str(output_dir / fname))
            # Return cached result with new output_dir
            cached["output_dir"] = str(output_dir)
            cached["cached"] = True
            cached["parsed_request"] = parsed_request.model_dump()
            cached["resolved_request"] = resolved_request.model_dump()
            # Update viz_html path to new output dir
            if cached.get("viz_html"):
                viz_new = output_dir / "route_view.html"
                cached["viz_html"] = str(viz_new) if viz_new.exists() else None
            (output_dir / "agent_result.json").write_text(
                json.dumps(cached, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            return cached

        # ── Run engine (cache miss) ──
        preferred = resolved_request.preferred_reactants
        if resolved_request.step_count == 1:
            planning_result = run_retrosynthesis(
                smiles=resolved_request.target_smiles,
                database_name=resolved_request.database,
                weights=resolved_request.weights,
                preferred_reactants=preferred,
                top_k=resolved_request.top_k,
                base_dir=str(ENGINE_ROOT),
                show_progress=False,
            )
        elif preferred:
            # Scenario C: multi-step with preferred reactant
            # Direct beam search with default 5 steps (no incremental trial)
            target_steps = resolved_request.step_count if resolved_request.step_count > 1 else 5
            planning_result = plan_retrosynthesis(
                smiles=resolved_request.target_smiles,
                steps=target_steps,
                database_name=resolved_request.database,
                weights=resolved_request.weights,
                preferred_reactants=preferred,
                base_dir=str(ENGINE_ROOT),
                beam_width=8,
                per_step_top_k=3,
            )
        else:
            planning_result = plan_retrosynthesis(
                smiles=resolved_request.target_smiles,
                steps=resolved_request.step_count,
                database_name=resolved_request.database,
                weights=resolved_request.weights,
                preferred_reactants=preferred,
                base_dir=str(ENGINE_ROOT),
                beam_width=8,
                per_step_top_k=3,
            )
        # ── Assign reaction types (cache lookup + LLM fallback) ──
        self._assign_reaction_types(planning_result)

        # ── Scenario C LLM fallback: if no engine route reached preferred reactants ──
        if preferred and self._scenario_c_needs_llm_fallback(planning_result):
            llm_routes = self._llm_design_routes(
                target_smiles=resolved_request.target_smiles,
                target_name=parsed_request.target_name or "",
                preferred_smiles_list=preferred,
                preferred_names=parsed_request.preferred_reactant_names,
                step_count=resolved_request.step_count,
                language=resolved_request.language,
            )
            if llm_routes:
                data = planning_result.setdefault("data", {})
                existing = data.get("all_routes", [])
                # Adjust ranks: engine routes first, LLM routes after
                next_rank = max((r.get("route_rank", 0) for r in existing), default=0) + 1
                for r in llm_routes:
                    r["route_rank"] = next_rank
                    next_rank += 1
                data["all_routes"] = existing + llm_routes
                # Update message
                planning_result["message"] = (
                    f"Engine found {len(existing)} route(s); "
                    f"LLM designed {len(llm_routes)} additional route(s)."
                )

        (output_dir / "planning_result.json").write_text(
            json.dumps(planning_result, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        # Generate visualization HTML
        viz_html_path = output_dir / "route_view.html"
        try:
            subprocess.run(
                [sys.executable, str(PROJECT_ROOT / "visualize.py"),
                 str(output_dir / "planning_result.json"),
                 "-o", str(viz_html_path),
                 "--lang", resolved_request.language],
                check=True, capture_output=True,
                cwd=str(PROJECT_ROOT),
            )
        except Exception:
            viz_html_path = None

        explanation = self.explain_results(user_message, resolved_request, planning_result)
        (output_dir / "agent_reply.md").write_text(explanation, encoding="utf-8")

        result = {
            "status": "success",
            "message": "Agent run completed.",
            "reply_markdown": explanation,
            "parsed_request": parsed_request.model_dump(),
            "resolved_request": resolved_request.model_dump(),
            "planning_result": planning_result,
            "output_dir": str(output_dir),
            "viz_html": str(viz_html_path) if viz_html_path and viz_html_path.exists() else None,
            "cached": False,
        }
        (output_dir / "agent_result.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        # ── Save to cache ──
        self._save_cache(cache_key, output_dir)

        return result

    # ── Cache helpers ──

    @staticmethod
    def _cache_key(resolved_request: ResolvedAgentRequest) -> str:
        """Build a deterministic cache key from target SMILES + scenario params."""
        from hashlib import md5
        smiles = resolved_request.target_smiles
        step = resolved_request.step_count
        db = resolved_request.database or "default"
        preferred = sorted(resolved_request.preferred_reactants)
        raw = f"{smiles}|{step}|{db}|{','.join(preferred)}"
        digest = md5(raw.encode()).hexdigest()[:8]
        # Use SMILES as readable prefix (sanitize for filesystem)
        safe_smiles = re.sub(r'[\\/:*?"<>|\s]', '_', smiles)[:40]
        return f"{safe_smiles}_{step}step_{db}_{digest}"

    @staticmethod
    def _load_cache(cache_key: str) -> Optional[Dict[str, Any]]:
        """Load cached result if available."""
        cache_dir = CACHE_ROOT / cache_key
        result_file = cache_dir / "agent_result.json"
        if not result_file.exists():
            return None
        try:
            data = json.loads(result_file.read_text(encoding="utf-8"))
            if data.get("status") == "success":
                return data
        except Exception:
            pass
        return None

    @staticmethod
    def _save_cache(cache_key: str, output_dir: Path) -> None:
        """Copy key output files to cache directory."""
        cache_dir = CACHE_ROOT / cache_key
        cache_dir.mkdir(parents=True, exist_ok=True)
        for fname in ("planning_result.json", "route_view.html", "agent_reply.md", "agent_result.json"):
            src = output_dir / fname
            if src.exists():
                shutil.copy2(str(src), str(cache_dir / fname))
