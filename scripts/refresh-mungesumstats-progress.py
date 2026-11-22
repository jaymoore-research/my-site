"""Publish fresh benchmark measurements to the GitHub Pages data file."""
import argparse
import base64
import json
import shlex
import subprocess
from pathlib import Path


def run(args, data=None):
    return subprocess.run(args, input=data, stdout=subprocess.PIPE, check=True, timeout=90).stdout


def add_github_availability(report):
    points = [point for panel in report["panels"] for point in panel["points"]]
    commits = sorted({point["commitSha"] for point in points})
    found = set()
    for start in range(0, len(commits), 25):
        block = commits[start:start + 25]
        fields = " ".join(f'c{i}: object(oid:"{sha}") {{ oid }}' for i, sha in enumerate(block))
        query = 'query { repository(owner:"neurogenomics", name:"mungesumstats-autoresearch") { ' + fields + ' } }'
        response = json.loads(run(["gh", "api", "graphql", "--input", "-"], json.dumps({"query": query}).encode()))
        if response.get("errors") or response["data"]["repository"] is None:
            raise SystemExit("Could not verify GitHub commit links")
        objects = response["data"]["repository"]
        found.update(sha for i, sha in enumerate(block) if (objects[f"c{i}"] or {}).get("oid") == sha)
    for point in points:
        point["github"] = point["commitSha"] in found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("host", help="SSH host where the benchmark receipts are stored")
    parser.add_argument("root", help="Absolute path to the benchmark directory on that host")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Compare data without publishing")
    mode.add_argument("--snapshot", type=Path, help="Save refreshed data locally without publishing")
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

    add_github_availability(report)
    fresh = (json.dumps(report, separators=(",", ":")) + "\n").encode()
    if args.snapshot:
        args.snapshot.write_bytes(fresh)
        return

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
