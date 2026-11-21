"""Publish fresh benchmark measurements to the GitHub Pages data file."""
import argparse
import base64
import json
import shlex
import subprocess
from pathlib import Path


def run(args, data=None):
    return subprocess.run(args, input=data, stdout=subprocess.PIPE, check=True, timeout=90).stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("host", help="SSH host where the benchmark receipts are stored")
    parser.add_argument("root", help="Absolute path to the benchmark directory on that host")
    parser.add_argument("--check", action="store_true", help="Compare data without publishing")
    args = parser.parse_args()
    exporter = Path(__file__).with_name("export-mungesumstats-progress.py")
    fresh = run(
        ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "--", args.host,
         "python3 - " + shlex.quote(args.root)],
        exporter.read_bytes(),
    )
    report = json.loads(fresh)
    if not report["panels"] or any(not panel["points"] for panel in report["panels"]):
        raise SystemExit("Benchmark export did not contain complete panels")

    endpoint = "repos/jaymoore-research/my-site/contents/mungesumstats/progress.json"
    live = json.loads(run(["gh", "api", endpoint]))
    if base64.b64decode(live["content"]) == fresh:
        print("Benchmark data is current.")
        return

    if args.check:
        print("New benchmark data is available.")
        return

    body = json.dumps({
        "message": "Refresh MungeSumstats benchmark measurements",
        "branch": "main",
        "sha": live["sha"],
        "content": base64.b64encode(fresh).decode("ascii"),
    }).encode("utf-8")
    receipt = json.loads(run(["gh", "api", "--method", "PUT", endpoint, "--input", "-"], body))
    print(f"Published benchmark data in {receipt['commit']['sha'][:12]}.")


if __name__ == "__main__":
    main()
