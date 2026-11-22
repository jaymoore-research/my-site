"""Read completed benchmark receipts and export public chart data."""
import csv
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from hashlib import sha256
from math import isfinite
from pathlib import Path
from statistics import median

SHA = re.compile(r"[0-9a-f]{40}")
TRIAL = re.compile(r"\d{8}T\d{6}-[0-9a-f]{8}")


def build(record, arm, stock):
    prefix = "baseline" if arm == "B" else "candidate"
    commit = record[f"{prefix}_commit"]
    dependency = record.get(f"{prefix}_deps", {}).get("BSgenome", stock)
    if not SHA.fullmatch(commit) or (dependency and not SHA.fullmatch(dependency)):
        raise ValueError("Invalid benchmark build ID")
    return commit, dependency


def measurements(studies, arm):
    runtime, ram = [], []
    for study in studies.values():
        runs = [run for run in study["runs"] if run["label"] == arm]
        if len(runs) != 2:
            raise ValueError("Expected two timed runs per build")
        seconds = [run["seconds"] for run in runs]
        rss = [run["peak_rss_kb"] for run in runs]
        if any(not isfinite(x) or x <= 0 for x in seconds + rss):
            raise ValueError("Invalid runtime or peak RSS in benchmark")
        runtime.append(median(seconds))
        ram.append(median(rss) / 1024**2)
    return round(median(runtime), 3), round(median(ram), 3)


def runtime_saved(studies):
    """Median of baseline minus candidate runtime within each input file, in seconds."""
    savings = []
    for study in studies.values():
        baseline, candidate = (
            median(run["seconds"] for run in study["runs"] if run["label"] == arm)
            for arm in ("B", "C")
        )
        savings.append(baseline - candidate)
    return round(median(savings), 3)


def export(root):
    cohorts = defaultdict(dict)
    baselines = {}
    for location in (root, root / "screen"):
        config = json.loads((location / "eval.json").read_text())
        stock = config.get("dependencies", {}).get("BSgenome", {}).get("stock_commit", "")
        with (location / "ledger.tsv").open() as source:
            rows = list(csv.DictReader(source, delimiter="\t"))
        promotions = {r["trial_id"] for r in rows if r.get("note") == "promoted"}
        for row in rows:
            trial = row.get("trial_id", "")
            if not TRIAL.fullmatch(trial) or row.get("gate") != "pass" or row["note"] == "promoted":
                continue
            receipt = location / "trials" / trial
            record = json.loads((receipt / "trial.json").read_text())
            for tier in ("2", "3", "holdout"):
                studies = record["tiers"].get(tier, {}).get("studies")
                if tier == "holdout" and isinstance(studies, int):
                    studies = json.loads((receipt / "holdout.json").read_text())["studies"]
                if not isinstance(studies, dict) or not studies:
                    continue
                kind = {"2": "Screening", "3": "Extended panel", "holdout": "Holdout"}[tier]
                parallel = all(s.get("side_by_side", False) for s in studies.values())
                names = tuple(sorted(studies))
                # Match the harness: reject a trial when more than one third of its studies are noisy.
                if sum(s["noisy"] for s in studies.values()) * 3 > len(studies):
                    continue
                try:
                    observed = {arm: measurements(studies, arm) for arm in ("B", "C")}
                except ValueError:
                    continue
                panel = kind, names, parallel
                group = cohorts[panel]
                status = "kept" if trial in promotions and row["tier_reached"] == tier else "tested"
                date = datetime.fromisoformat(row["at_utc"]).isoformat(timespec="seconds")
                baseline = build(record, "B", stock)
                if panel not in baselines or date < baselines[panel]["date"]:
                    commit, _ = baseline
                    seconds, ram = observed["B"]
                    baselines[panel] = dict(commit=commit[:12], date=date, seconds=seconds, gib=ram,
                                            stockDependency=stock[:12])
                for arm in ("B", "C"):
                    key = build(record, arm, stock)
                    point_status = status if arm == "C" else "reference"
                    old = group.get(key)
                    if old and point_status == "reference":
                        continue
                    if old and ({"reference":0,"tested":1,"kept":2}[old["status"]], old["date"]) >= (
                        {"reference":0,"tested":1,"kept":2}[point_status], date
                    ):
                        continue
                    runtime, ram = observed[arm]
                    change = None
                    if arm == "C":
                        mss_changed, bsg_changed = (a != b for a, b in zip(key, baseline))
                        change = "both" if mss_changed and bsg_changed else (
                            "mss" if mss_changed else "bsg" if bsg_changed else "same"
                        )
                    group[key] = dict(
                        commit=key[0][:12], dependency=key[1][:12], date=date,
                        commitSha=key[0], dependencySha=key[1], change=change,
                        pairedBase=baseline[0][:12] if change else None,
                        savedSeconds=runtime_saved(studies) if change else None,
                        first=min(date, old["first"]) if old else date, measured=arm,
                        status=point_status, seconds=runtime, gib=ram,
                    )
    panels = []
    for (kind, names, parallel), group in cohorts.items():
        points = sorted(group.values(), key=lambda p: p["first"])
        # Initial baseline and candidate can share a timestamp. Insertion order puts the baseline first.
        identity = sha256(repr((kind, names, parallel)).encode()).hexdigest()[:10]
        panels.append(dict(id=identity, kind=kind, studies=len(names), parallel=parallel,
                           baseline=baselines[kind, names, parallel], points=points))
    order = {"Screening":0,"Holdout":1,"Extended panel":2}
    panels.sort(key=lambda p: (order[p["kind"]], -len(p["points"]), p["id"]))
    if not panels:
        raise SystemExit("No completed timing receipts found")
    return dict(panels=panels)


if __name__ == "__main__":
    json.dump(export(Path(sys.argv[1])), sys.stdout, separators=(",", ":"))
    print()
