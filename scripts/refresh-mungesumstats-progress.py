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


def publish(path, fresh, message, check):
    """Commit one data file to my-site main through the contents API, unless it is current."""
    endpoint = f"repos/jaymoore-research/my-site/contents/{path}"
    try:
        live = json.loads(run(["gh", "api", endpoint]))
        current, sha = base64.b64decode(live["content"]), live["sha"]
    except subprocess.CalledProcessError:
        current, sha = None, None  # a new file
    if current == fresh:
        print(f"{path} is current.")
        return
    if check:
        print(f"New data for {path} is available.")
        return
    body = {"message": message, "branch": "main", "content": base64.b64encode(fresh).decode("ascii")}
    if sha: body["sha"] = sha
    receipt = json.loads(run(["gh", "api", "--method", "PUT", endpoint, "--input", "-"], json.dumps(body).encode()))
    print(f"Published {path} in {receipt['commit']['sha'][:12]}.")


def export_dbsnp(host, root):
    """The dbSNP lane's tab data. A missing lane is not an error for the loop's data."""
    exporter = Path(__file__).with_name("export-dbsnp-progress.py")
    report = json.loads(run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "--", host,
                             "python3 - " + shlex.quote(root)], exporter.read_bytes()))
    # The generated time changes every run; leave it out of the comparison by keeping it only when data changed.
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("host", help="SSH host where the benchmark receipts are stored")
    parser.add_argument("root", help="Absolute path to the benchmark directory on that host")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Compare data without publishing")
    mode.add_argument("--snapshot", type=Path, help="Save refreshed data locally without publishing")
    parser.add_argument("--dbsnp-root", default="/mnt/data3/jmoore3/mss-dbsnp",
                        help="dbSNP lane root on the same host; empty string skips the dbSNP tab")
    args = parser.parse_args()
    if args.dbsnp_root:
        lane = export_dbsnp(args.host, args.dbsnp_root)
        endpoint = "repos/jaymoore-research/my-site/contents/mungesumstats/dbsnp.json"
        try:
            old = json.loads(base64.b64decode(json.loads(run(["gh", "api", endpoint]))["content"]))
            if {k: v for k, v in old.items() if k != "generated"} == {k: v for k, v in lane.items() if k != "generated"}:
                lane["generated"] = old["generated"]
        except (subprocess.CalledProcessError, ValueError, KeyError):
            pass
        publish("mungesumstats/dbsnp.json", (json.dumps(lane, separators=(",", ":")) + "\n").encode(),
                "Refresh dbSNP lane progress", args.check or bool(args.snapshot))
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
