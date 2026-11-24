"""Read the dbSNP lane's ledger and state and export public data for the dbSNP tab.

    python3 - /mnt/data3/jmoore3/mss-dbsnp [/mnt/data3/jmoore3/mss-cx3/mss-dbsnp] < export-dbsnp-progress.py

Each root is one machine's copy of the lane: iclneuro's own, and the CX3 results that pull-results.sh
copies to iclneuro. Turns from all roots are merged; branch state and the baseline come from the first
root that has them, so list iclneuro first: turns are timed there.

The lane writes holdout inputs as holdout-N, so no accession reaches this file. Rows whose detail
still names a GCST accession are refused rather than published.
"""
import csv
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

# Branch, the backend directory it measures, and a plain description. Order is the page's order.
BRANCHES = [
    ("snplocs-baseline", "snplocs", "SNPlocs packages, as MungeSumstats uses them now (comparator)"),
    ("build", "build", "Release builder: NCBI files to one shared table"),
    ("duckdb-parquet", "duckdb-parquet", "Parquet sorted by id and by position, queried with DuckDB"),
    ("duckdb-native", "duckdb-native", "One DuckDB file per genome build"),
    ("arrow-parquet", "arrow-parquet", "Parquet through R arrow (tidyGWAS layout)"),
    ("flat-mmap", "flat-mmap", "Sorted binary files merged in C over a memory map"),
    ("hdf5", "hdf5", "HDF5 through rhdf5 (gwaslab layout)"),
    ("sqlite", "sqlite", "SQLite tables keyed by id and by position"),
    ("genomicscores", "genomicscores", "GenomicScores layout: per-position tracks, exact one-byte allele code"),
    ("vortex", "vortex", "Vortex files through DuckDB"),
    ("hail", "hail", "Hail tables through BiocHail"),
    ("echtvar-l2", "echtvar-l2", "Position-and-allele lookup with echtvar"),
]
ACCESSION = re.compile(r"GCST\d+", re.I)
NUMBER = ("cold_s", "cold_gb", "warm_s", "warm_gb", "ref_cold_s", "ref_cold_gb", "time_ratio", "ram_ratio")


def number(value):
    try:
        return round(float(value), 4)
    except (TypeError, ValueError):
        return None


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def machine(root):
    return "CX3" if "mss-cx3" in str(root) else "iclneuro"


def first(roots, rel):
    for root in roots:
        found = read_json(root / rel)
        if found: return found, machine(root)
    return None, None


def export(*roots):
    rows = []
    for root in roots:
        ledger = root / "records" / "dbsnp-ledger.tsv"
        if ledger.exists():
            rows += [dict(r, machine=machine(root)) for r in csv.DictReader(ledger.open(), delimiter="\t")]
    rows.sort(key=lambda r: r.get("utc") or "")
    # Holdout accessions, read on the server so they can be refused, never written out.
    if rows and not any(p for root in roots for p in (root / "hidden").glob("*.rds")):
        raise SystemExit("No holdout list beside the ledger, so holdout names cannot be refused; refusing to publish")
    hidden = {m.group(0).upper() for root in roots for p in (root / "hidden").glob("*.rds") for m in ACCESSION.finditer(p.stem)}
    turns = []
    for row in rows:
        detail = row.get("b0_detail", "")
        text = " ".join(str(v) for v in row.values()).upper()
        if any(acc in text for acc in hidden):
            raise SystemExit("A ledger row names a holdout study; refusing to publish")
        if re.search(r"/mnt/|/rds/", detail):
            raise SystemExit("Ledger detail names a private path; refusing to publish")
        turn = {k: row.get(k) or None for k in ("utc", "backend", "commit", "kind", "release", "b0", "ref", "verdict", "note", "machine")}
        turn.update({k: number(row.get(k)) for k in NUMBER})
        turn["b0Detail"] = detail[:300]
        turns.append(turn)

    base, base_machine = first(roots, "state/snplocs/155.screen.json")
    baseline = None
    if base:
        baseline = dict(cold_s=round(median(base["cold_s"]), 3), cold_gb=round(median(base["cold_gb"]), 3),
                        warm_s=round(median(base["warm_s"]), 3), warm_gb=round(median(base["warm_gb"]), 3), machine=base_machine)

    branches = []
    for branch, backend, description in BRANCHES:
        mine = [t for t in turns if t["backend"] == backend]
        screens = [t for t in mine if t["kind"] == "screen"]
        gates = [t for t in mine if t["kind"] == "gate-a"]
        best, _ = first(roots, f"state/{backend}/155.best.json")
        entry, _ = first(roots, f"state/{backend}/entry.json")
        if gates:
            stage = "Gate A " + (gates[-1]["verdict"] or "")
        elif screens:
            stage = "Turns"
        elif entry and entry.get("open"):
            stage = "Ready"
        else:
            stage = "Waiting"
        branches.append(dict(
            branch="dbsnp/" + branch, backend=backend, description=description, stage=stage,
            entryOpen=bool(entry and entry.get("open")),
            entryFailed=sorted({f["id"] for f in (entry or {}).get("failed", [])}),
            turns=len(screens), better=sum(t["verdict"] == "better" for t in screens),
            lastVerdict=screens[-1]["verdict"] if screens else None,
            bestColdSeconds=round(median(best["cold_s"]), 3) if best else None,
            bestColdGB=round(median(best["cold_gb"]), 3) if best else None,
            gateA=gates[-1]["verdict"] if gates else None,
        ))
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return dict(generated=generated, release="155", input="1M random rs ids, GRCh38, one core",
                baseline=baseline, branches=branches, turns=turns)


if __name__ == "__main__":
    json.dump(export(*(Path(a) for a in sys.argv[1:])), sys.stdout, separators=(",", ":"))
    print()
