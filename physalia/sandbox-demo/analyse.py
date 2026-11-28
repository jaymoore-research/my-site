"""Summarise group measurements from samples.json into summary.json."""
import json
from pathlib import Path

root = Path(__file__).resolve().parent
rows = json.loads((root / "samples.json").read_text())
groups = {}
for row in rows:
    groups.setdefault(row["group"], []).append(row["value"] or 0)
summary = {name: {"mean": sum(values) / len(values), "n": len(values)}
           for name, values in groups.items()}
(root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print("Wrote summary.json.")
