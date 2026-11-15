"""
Build manifest
==============

What this does
--------------
Scans the repository for data files and writes manifest.json, a short index of every configuration that the
viewer can read instead of asking GitHub for the file list. This step is optional. Without a manifest, the
viewer finds the data files on its own.

How to run
----------
Run it from the repository root after adding or moving data files.

    python3 build_manifest.py            write manifest.json
    python3 build_manifest.py --strip    also shrink finished data files

With --strip, finished data files lose their "_checked" list, which is the closure search's record of settled
pairs and about 75% of a large file. They are tagged "_complete": true instead, so the search knows not to resume
them. Unfinished files and compressed (.json.gz) files are left untouched.

For reviewers
-------------
Data files are data_<configuration>.json or data_<configuration>.json.gz. Rank files (rank_*.json) and exact
intersection reports (closure_*.json) are attached to their configuration when present.
"""

import os
import sys
import gzip
import json

ROOT = os.path.dirname(os.path.abspath(__file__))


def row_lengths(config):
    """The row lengths of a configuration. For example, '332' -> [3, 3, 2] and '1-10' -> [1, 10]."""
    if "-" in config:
        return [int(x) for x in config.split("-")]
    return [int(d) for d in config]


def walk_files():
    """Every file in the repository as (folder, file name), skipping hidden folders such as .git."""
    for folder, subfolders, files in os.walk(ROOT):
        subfolders[:] = sorted(d for d in subfolders if not d.startswith("."))
        for name in sorted(files):
            yield folder, name


def relative(folder, name):
    return os.path.relpath(os.path.join(folder, name), ROOT).replace(os.sep, "/")


def config_of(name, prefix):
    """'data_332.json.gz' -> '332', or None if the name does not match the prefix."""
    for suffix in (".json.gz", ".json"):
        if name.startswith(prefix) and name.endswith(suffix):
            return name[len(prefix) : -len(suffix)]
    return None


def load(path):
    opener = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with opener as f:
        return json.load(f)


def strip_checked(path, data):
    """Remove the list of settled pairs from a finished file and return how many bytes this saved."""
    before = os.path.getsize(path)
    data.pop("_checked")
    data["_complete"] = True
    temporary = path + ".tmp"
    with open(temporary, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(temporary, path)
    return before - os.path.getsize(path)


def main():
    strip = "--strip" in sys.argv

    # First pass: rank files and intersection reports, so each one can be attached to its configuration.
    reports, ranks = {}, {}
    for folder, name in walk_files():
        if config_of(name, "closure_"):
            reports[config_of(name, "closure_")] = relative(folder, name)
        if config_of(name, "rank_"):
            ranks[config_of(name, "rank_")] = relative(folder, name)

    # Second pass: one entry per data file. A compressed file wins over a plain one for the same configuration.
    entries, stripped, saved = {}, 0, 0
    for folder, name in walk_files():
        if not config_of(name, "data_"):
            continue
        path = os.path.join(folder, name)
        data = load(path)
        config = data["partition"]
        if config in entries and not name.endswith(".gz"):
            continue
        rows = row_lengths(config)
        num_cells = len(data["cells"])
        complete = data.get("_complete") or len(data.get("_checked", [])) >= num_cells * (num_cells - 1)
        if strip and complete and "_checked" in data and not name.endswith(".gz"):
            saved += strip_checked(path, data)
            stripped += 1
        dims = [cell["dim"] for cell in data["cells"]]
        entry = {
            "partition": config,
            "rows": rows,
            "n": sum(rows),
            "path": relative(folder, name),
            "cells": num_cells,
            "edges": len(data["edges"]),
            "maxDim": max(dims),
            "dimProfile": [dims.count(k) for k in range(max(dims) + 1)],
            "labeling": data.get("labeling"),
            "complete": bool(complete),
            "bytes": os.path.getsize(path),
        }
        if config in reports:
            entry["report"] = reports[config]
        if config in ranks:
            entry["rank"] = ranks[config]
        entries[config] = entry

    ordered = sorted(entries.values(), key=lambda e: (len(e["rows"]), e["n"], e["rows"]))
    with open(os.path.join(ROOT, "manifest.json"), "w") as f:
        json.dump({"partitions": ordered}, f, indent=1)
    print(f"Wrote manifest.json with {len(ordered)} configurations.")
    if strip:
        print(f"Shrank {stripped} finished data files, saving {saved / 1e6:.1f} MB.")


if __name__ == "__main__":
    main()
