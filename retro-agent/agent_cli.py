import argparse
import json
import sys

from agent_runtime import RetrosynthesisAgent

# Fix Unicode rendering on Windows consoles (GBK codec can't encode emoji / CJK ext)
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Natural-language retrosynthesis agent CLI")
    parser.add_argument("-q", "--query", type=str, help="Natural-language user request")
    parser.add_argument("--model", type=str, default=None, help="Override the OpenAI model name")
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the full JSON result instead of only the markdown reply",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    query = args.query or input("User request: ").strip()
    if not query:
        print("No user request provided.")
        return 1

    agent = RetrosynthesisAgent(model=args.model)
    result = agent.run(query)
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print(result["reply_markdown"])
        print(f"\n[output_dir] {result['output_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
