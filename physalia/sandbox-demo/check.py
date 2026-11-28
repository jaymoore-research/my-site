"""Check the expected output for this fixed synthetic teaching fixture."""
import json
import sys
from pathlib import Path

expected = {"control": {"mean": 12.0, "n": 2},
            "treated": {"mean": 14.0, "n": 3}}
try:
    actual = json.loads((Path(__file__).resolve().parent / "summary.json").read_text())
except (OSError, ValueError) as error:
    print(f"FAIL: summary.json is missing or unreadable: {error}")
    sys.exit(1)
if actual != expected:
    print(f"FAIL: expected {expected}; found {actual}")
    sys.exit(1)
print("PASS: group means and measured sample counts match the teaching fixture.")
