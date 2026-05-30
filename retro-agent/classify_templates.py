#!/usr/bin/env python3
"""
One-time script: classify all SMARTS reaction templates using LLM.

Reads SimpRetro4Learn/reaction_template.json, sends unique templates in batches
to an LLM for reaction type classification, and saves the mapping to
template_reaction_types.json for fast runtime lookup.

Usage:
  python classify_templates.py                      # classify all, default batch_size=25
  python classify_templates.py --batch-size 30      # custom batch size
  python classify_templates.py --resume             # skip already-classified templates
  python classify_templates.py --dry-run            # show prompt for first batch, no API call

Requires: OPENAI_API_KEY and optionally OPENAI_BASE_URL in environment.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

# ── Paths ──
AGENT_ROOT = Path(__file__).resolve().parent
ENGINE_ROOT = AGENT_ROOT.parent / "SimpRetro4Learn"
TEMPLATE_PATH = ENGINE_ROOT / "reaction_template.json"
OUTPUT_PATH = AGENT_ROOT / "template_reaction_types.json"


def load_unique_templates():
    """Load all unique SMARTS templates from reaction_template.json."""
    if not TEMPLATE_PATH.exists():
        print(f"Error: {TEMPLATE_PATH} not found")
        sys.exit(1)
    with open(TEMPLATE_PATH, encoding="utf-8") as f:
        raw = json.load(f)
    # Deduplicate while preserving order
    seen = set()
    unique = []
    for t in raw:
        t = t.strip()
        if t and t not in seen:
            seen.add(t)
            unique.append(t)
    return unique


def load_existing(output_path):
    """Load previously classified templates for resume support."""
    if not output_path.exists():
        return {}
    with open(output_path, encoding="utf-8") as f:
        return json.load(f)


def setup_llm_client():
    """Create an OpenAI-compatible client from environment variables."""
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("Error: OPENAI_API_KEY is not set in environment.")
        print("Set it before running: $env:OPENAI_API_KEY='your-key'")
        sys.exit(1)

    from openai import OpenAI
    base_url = os.getenv("OPENAI_BASE_URL")
    if base_url:
        return OpenAI(api_key=api_key, base_url=base_url)
    return OpenAI(api_key=api_key)


def classify_batch(client, templates, batch_start, model=None, dry_run=False):
    """Send one batch of templates to LLM for classification.

    Returns list of {"template_index": N, "reaction_type": "..."} dicts.
    """
    batch_items = []
    for i, t in enumerate(templates):
        idx = batch_start + i
        batch_items.append(f"T{idx}: {t}")

    prompt = f"""You are an expert organic chemist. Classify each SMARTS reaction template below.

Each template is a SMARTS string in the format REACTANTS>>PRODUCTS:
- >> separates reactants (left side, what is consumed) from products (right side, what is formed)
- . separates multiple molecules on the same side
- [C:1], [O:2] etc. are atom-mapped atoms showing how atoms move
- [C;H0;D3;+0:1] means carbon with 0 H, degree 3, charge 0, atom map 1

Read the SMARTS to understand what transformation each template encodes, then classify.

Classification categories (choose the MOST specific applicable):
- "Diels-Alder cycloaddition" — two alkenes/conjugated system → cyclohexene ring
- "Cross-coupling" — organometallic + organic halide → new C-C bond
- "Esterification" — acid + alcohol → ester
- "Hydrolysis" — ester/amide/nitrile + H2O → acid + alcohol/amine
- "Oxidation" — gain of O, loss of H (e.g. alcohol→carbonyl, alkene→epoxide)
- "Reduction" — gain of H, loss of O (e.g. carbonyl→alcohol)
- "Elimination" — formation of C=C with loss of H2O/HX/halogen
- "Substitution (SN2)" — bimolecular nucleophilic substitution
- "Substitution (SN1)" — unimolecular nucleophilic substitution
- "Grignard addition" — RMgX + carbonyl → alcohol
- "Wittig reaction" — phosphonium ylide + carbonyl → alkene
- "Aldol condensation" — enolate + carbonyl → alpha,beta-unsaturated carbonyl
- "Catalytic hydrogenation" — H2 addition to C=C/C≡C (often has Pd/Pt/Ni in template or conditions)
- "Halogenation" — addition of Br2/Cl2 across C=C, or radical halogenation
- "Alkylation" — introduction of alkyl group (e.g. Williamson ether synthesis)
- "Epoxidation" — alkene → epoxide
- "Amidation" — acid/ester + amine → amide
- "Protection" — introduction of protecting group
- "Deprotection" — removal of protecting group
- "Rearrangement" — skeletal rearrangement without adding/removing atoms
- "Other" — none of the above clearly apply

Return JSON only:
{{"classifications": [{{"template_index": {batch_start}, "reaction_type": "Diels-Alder cycloaddition"}}, ...]}}

Templates to classify:
{chr(10).join(batch_items)}"""

    if dry_run:
        print(f"\n{'='*60}")
        print(f"DRY RUN — Batch {batch_start // len(templates) + 1}")
        print(f"Templates {batch_start}–{batch_start + len(templates) - 1}")
        print(f"{'='*60}")
        print(f"Prompt length: {len(prompt)} chars")
        print(f"\nFirst 2000 chars of prompt:")
        print(prompt[:2000])
        return []

    model_name = model or os.getenv("OPENAI_MODEL", "gpt-5.4-mini")
    # Check if using DeepSeek API
    base_url = os.getenv("OPENAI_BASE_URL", "").lower()
    is_deepseek = "deepseek" in model_name.lower() or "deepseek.com" in base_url

    if is_deepseek:
        response = client.chat.completions.create(
            model=model_name,
            messages=[
                {"role": "system", "content": "Return valid json only."},
                {"role": "user", "content": prompt},
            ],
            response_format={"type": "json_object"},
            stream=False,
        )
        content = response.choices[0].message.content or "{}"
    else:
        response = client.responses.create(
            model=model_name,
            input=prompt,
        )
        content = response.output_text

    # Parse JSON
    content = content.strip()
    if content.startswith("```"):
        import re
        content = re.sub(r"^```(?:json)?", "", content).strip()
        content = re.sub(r"```$", "", content).strip()
    try:
        result = json.loads(content)
        return result.get("classifications", [])
    except json.JSONDecodeError:
        print(f"  Warning: could not parse JSON response for batch starting at {batch_start}")
        print(f"  Response (first 500 chars): {content[:500]}")
        return []


def main():
    parser = argparse.ArgumentParser(description="Classify all SMARTS reaction templates via LLM")
    parser.add_argument("--batch-size", type=int, default=25, help="Templates per LLM call (default: 25)")
    parser.add_argument("--resume", action="store_true", help="Skip already-classified templates")
    parser.add_argument("--dry-run", action="store_true", help="Show prompt for first batch only, no API call")
    parser.add_argument("--model", type=str, default=None, help="Model name (default: OPENAI_MODEL env or gpt-5.4-mini)")
    parser.add_argument("--delay", type=float, default=1.5, help="Seconds between batches (default: 1.5)")
    args = parser.parse_args()

    # Load templates
    all_templates = load_unique_templates()
    print(f"Loaded {len(all_templates)} unique templates from {TEMPLATE_PATH}")

    # Resume: load existing classifications
    existing = {}
    if args.resume and OUTPUT_PATH.exists():
        existing = load_existing(OUTPUT_PATH)
        print(f"Resuming: {len(existing)} templates already classified")

    # Filter out already-classified
    to_classify = []
    for i, t in enumerate(all_templates):
        if t not in existing:
            to_classify.append((i, t))

    if not to_classify:
        print("All templates already classified. Nothing to do.")
        return

    print(f"Need to classify: {len(to_classify)} templates "
          f"({len(existing)} already done, {len(all_templates) - len(existing) - len(to_classify)} skipped)")

    if args.dry_run:
        # Show first batch prompt
        batch = [t for _, t in to_classify[:args.batch_size]]
        classify_batch(None, batch, 0, model=args.model, dry_run=True)
        return

    # Setup LLM client
    client = setup_llm_client()
    model_name = args.model or os.getenv("OPENAI_MODEL", "gpt-5.4-mini")
    print(f"Using model: {model_name}")

    # Classify in batches
    results = dict(existing)  # start with existing classifications
    total_batches = (len(to_classify) + args.batch_size - 1) // args.batch_size

    for batch_num in range(total_batches):
        start_idx = batch_num * args.batch_size
        end_idx = min(start_idx + args.batch_size, len(to_classify))
        batch_items = to_classify[start_idx:end_idx]
        batch_templates = [t for _, t in batch_items]

        print(f"\nBatch {batch_num + 1}/{total_batches}: "
              f"classifying templates {batch_items[0][0]}–{batch_items[-1][0]} "
              f"({len(batch_items)} templates)...")

        try:
            classifications = classify_batch(
                client, batch_templates,
                batch_start=batch_items[0][0],
                model=model_name,
            )
        except Exception as e:
            print(f"  Error: {e}")
            print(f"  Saving progress and exiting...")
            break

        # Store results
        for c in classifications:
            idx = c.get("template_index", -1)
            rxn_type = c.get("reaction_type", "")
            if isinstance(idx, int) and 0 <= idx < len(all_templates):
                results[all_templates[idx]] = rxn_type

        # Save after each batch (incremental)
        with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

        batch_classified = len(classifications)
        total_classified = len(results)
        print(f"  Classified {batch_classified} in this batch "
              f"({total_classified}/{len(all_templates)} total)")

        # Delay between batches to respect rate limits
        if batch_num < total_batches - 1:
            time.sleep(args.delay)

    print(f"\nDone! {len(results)}/{len(all_templates)} templates classified.")
    print(f"Saved to: {OUTPUT_PATH}")

    # Show some stats
    if results:
        type_counts = {}
        for rxn_type in results.values():
            type_counts[rxn_type] = type_counts.get(rxn_type, 0) + 1
        print("\nTop reaction types:")
        for rxn_type, count in sorted(type_counts.items(), key=lambda x: -x[1])[:15]:
            print(f"  {count:4d}  {rxn_type}")


if __name__ == "__main__":
    main()
