"""
Rank test
=========

What this does
--------------
The closure search finds a relation "the closure of C_w meets C_v" by exhibiting a curve. The rank test is a
quick check from the other side. Ranks can only drop in a limit, so

    the closure of C_w can meet C_v only if v <= w in Bruhat order, and
    r_k(i, j) = dim(X^k V_j + V_i) - i at the flag of v is at most its value at a generic point of C_w,
    for all i, j and k.

For each configuration, this script lists the pairs that the test allows, compares them with the relations the
search found, and writes the result to rank_<configuration>.json next to the data file. If the data is consistent,
no found relation is ever forbidden by the test.

How to run
----------
Run it from the repository root.

    python3 rank_test.py                 every configuration
    python3 rank_test.py 333             only (3,3,3)
    python3 rank_test.py 333 532 1-10    several configurations

Settings (environment variables, all optional)
----------------------------------------------
    REPO=<folder>      where the repository is (default: the current folder)
    OUT=<folder>       where to write the rank files (default: next to each data file)
    WORKERS=2          number of processes

For reviewers
-------------
All linear algebra is exact, modulo the prime P. The generic value of each rank is taken at one random point of
the cell. This gives the generic rank unless the point happens to land on a lower-dimensional subset, which has
negligible probability for a prime this large.
"""

import os
import sys
import json
import glob
import gzip
import time
import random
import multiprocessing as mp
from itertools import combinations

import numpy as np

P = 2147483647  # every computation is done modulo this prime

REPO = os.environ.get("REPO", ".")
OUT = os.environ.get("OUT", REPO)
WORKERS = int(os.environ.get("WORKERS", "2"))


# ── Reading the data ─────────────────────────────────────────────────────────────────────────────────────────────


def parse_perm(text):
    """'3621' -> [3, 6, 2, 1]. Here, '10' counts as one entry, and dashes separate entries once n >= 10."""
    if "-" in text:
        return [int(x) for x in text.split("-")]
    out, i = [], 0
    while i < len(text):
        if i + 1 < len(text) and text[i + 1] == "0":
            out.append(int(text[i : i + 2]))
            i += 2
        else:
            out.append(int(text[i]))
            i += 1
    return out


def rows_of(config):
    """The row lengths of a configuration. For example, '332' -> [3, 3, 2] and '1-10' -> [1, 10]."""
    return [int(x) for x in (config.split("-") if "-" in config else config)]


def tableau(rows):
    """The boxes labeled 1, ..., n down each column, starting from the leftmost column."""
    labels = [[0] * length for length in rows]
    label = 1
    for column in range(max(rows)):
        for row in labels:
            if column < len(row):
                row[column] = label
                label += 1
    return labels


def right_neighbor(rows, n):
    """right[b] is the box right of box b (counting from 0), or -1 at the end of a row. So (Xv)[b] = v[right[b]]."""
    right = [-1] * n
    for row in tableau(rows):
        for k in range(1, len(row)):
            right[row[k - 1] - 1] = row[k] - 1
    return right


COMPILED = {}  # each Plücker coordinate is compiled once and then reused


def evaluate(expression, values):
    """Evaluate a Plücker coordinate such as 'a3 - a1*a2' at the given values."""
    code = COMPILED.get(expression)
    if code is None:
        code = COMPILED[expression] = compile(expression, "<plucker>", "eval")
    return eval(code, {"__builtins__": {}}, values)


# ── Linear algebra modulo P ──────────────────────────────────────────────────────────────────────────────────────


def cell_matrix(cell, n, values, index):
    """
    Rebuild the matrix C of a cell at the given values of its free variables, from the stored Plücker coordinates.
    In other words, C[s][col] = ±p(pivots of the earlier columns, together with s).
    """
    w = parse_perm(cell["perm"])
    plucker = [level.strip("[]").split(" : ") for level in cell["plucker"]]
    C = [[0] * n for _ in range(n)]
    for col in range(n - 1):
        pivots = [w[c] - 1 for c in range(col)]
        for s in range(n):
            if s in pivots:
                continue
            rows = sorted(pivots + [s])
            position = rows.index(s)
            order = [w.index(p + 1) for p in rows if p != s]
            inversions = sum(1 for a in range(len(order)) for b in range(a + 1, len(order)) if order[a] > order[b])
            value = evaluate(plucker[col][index[col + 1][tuple(rows)]].strip(), values) % P
            C[s][col] = value if (position + col + inversions) % 2 == 0 else (-value) % P
    C[w[n - 1] - 1][n - 1] = 1
    return C, w


def rank_profile(C, w, n, right, K):
    """Every r_k(i, j), for k = 1, ..., K, i = 1, ..., n-1 and j = 1, ..., n, as one flat list."""
    columns = [[C[r][c] for r in range(n)] for c in range(n)]
    bottom_up = sorted(range(n), key=lambda c: -w[c])
    out = []
    Y = columns
    for k in range(1, K + 1):
        Y = [[y[right[r]] if right[r] >= 0 else 0 for r in range(n)] for y in Y]  # apply X once more
        for i in range(1, n):
            V_i = [c for c in bottom_up if c < i]  # the first i columns of C
            basis = {}  # pivot row -> reduced vector
            rank = 0
            for j in range(n):
                y = Y[j][:]
                for c in V_i:  # reduce modulo V_i
                    t = y[w[c] - 1]
                    if t:
                        y = [(a - t * b) % P for a, b in zip(y, columns[c])]
                for pivot, b in basis.items():  # reduce modulo the earlier columns of X^k V_j
                    t = y[pivot]
                    if t:
                        y = [(a - t * bb) % P for a, bb in zip(y, b)]
                lead = next((r for r in range(n) if y[r]), None)
                if lead is not None:
                    inverse = pow(y[lead], P - 2, P)
                    y = [(a * inverse) % P for a in y]
                    for pivot in list(basis):
                        t = basis[pivot][lead]
                        if t:
                            basis[pivot] = [(a - t * bb) % P for a, bb in zip(basis[pivot], y)]
                    basis[lead] = y
                    rank += 1
                out.append(rank)
    return out


def bruhat_matrix(w, n):
    """The Bruhat rank matrix, r[i][j] = #{a <= j : w(a) <= i}."""
    A = np.zeros((n, n), dtype=np.int16)
    for column, value in enumerate(w):
        A[value - 1, column] = 1
    return A.cumsum(0).cumsum(1)


# ── The test for one configuration ───────────────────────────────────────────────────────────────────────────────


def run(path, rng):
    started = time.time()
    data = json.load(gzip.open(path, "rt") if path.endswith(".gz") else open(path))
    config = data["partition"]
    rows = rows_of(config)
    n = sum(rows)
    # A finished search has settled every ordered pair of cells. Older finished files have no '_complete' tag,
    # so an unfinished file is one whose list of settled pairs is still shorter than that.
    num_cells = len(data["cells"])
    settled = len(set(data.get("_checked", [])))
    complete = bool(data.get("_complete")) or "_checked" not in data or settled >= num_cells * (num_cells - 1)
    if not data["edges"] and not data.get("_complete"):
        return config, None, complete

    right = right_neighbor(rows, n)
    K = max(rows) - 1
    index = [None] + [{rs: j for j, rs in enumerate(combinations(range(n), r))} for r in range(1, n)]
    cells = data["cells"]
    N = len(cells)

    generic = np.zeros((N, K * (n - 1) * n), dtype=np.int16)  # ranks at a random point of each cell
    at_flag = np.zeros_like(generic)  # ranks at the permutation flag, where every variable is 0
    bruhat = np.zeros((N, n, n), dtype=np.int16)
    for a, cell in enumerate(cells):
        C, w = cell_matrix(cell, n, {v: rng.randrange(1, P) for v in cell["free"]}, index)
        generic[a] = rank_profile(C, w, n, right, K)
        C0, _ = cell_matrix(cell, n, {v: 0 for v in cell["free"]}, index)
        at_flag[a] = rank_profile(C0, w, n, right, K)
        bruhat[a] = bruhat_matrix(w, n)
    bruhat = bruhat.reshape(N, -1)

    position_of = {cell["perm"]: i for i, cell in enumerate(cells)}
    found = [set() for _ in range(N)]
    for edge in data["edges"]:
        found[position_of[edge["src"]]].add(position_of[edge["tgt"]])

    extra, forbidden, allowed_total, closure_sizes = {}, 0, 0, []
    for s in range(N):
        allowed_mask = np.all(bruhat >= bruhat[s], axis=1) & np.all(at_flag <= generic[s], axis=1)
        allowed_mask[s] = False
        allowed = set(np.nonzero(allowed_mask)[0].tolist())
        allowed_total += len(allowed)
        forbidden += len(found[s] - allowed)
        no_curve_yet = sorted(allowed - found[s])
        if no_curve_yet:
            extra[cells[s]["perm"]] = [cells[t]["perm"] for t in no_curve_yet]
        closure_sizes.append(len(allowed) + 1)

    dims = [cell["dim"] for cell in cells]
    two_row = sum(closure_sizes[s] == 2 ** dims[s] for s in range(N)) if len(rows) == 2 else None
    summary = {
        "cells": N,
        "found": len(data["edges"]),
        "allowed": allowed_total,
        "predicted_only": sum(len(v) for v in extra.values()),
        "found_but_forbidden": forbidden,
        "two_row_closures_with_2_to_dim": two_row,
    }
    out = {
        "partition": config,
        "labeling": data.get("labeling"),
        "method": "rank test: v <= w in Bruhat order and dim(X^k V_j + V_i) at the "
        "permutation flag of v <= its generic value on C_w, for all i, j, k",
        "summary": summary,
        "extra": extra,
    }

    destination = os.path.join(OUT, os.path.dirname(os.path.relpath(path, REPO)), f"rank_{config}.json")
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    json.dump(out, open(destination, "w"), separators=(",", ":"))
    return config, dict(summary, secs=round(time.time() - started, 1)), complete


def job(path):
    seed = sum(map(ord, os.path.basename(path)))  # a fixed seed per file, so every run gives the same output
    return run(path, random.Random(seed))


def config_key(path):
    """'.../data_332.json.gz' -> '332'."""
    return os.path.basename(path)[len("data_") :].split(".json")[0]


# ── Printing the results ─────────────────────────────────────────────────────────────────────────────────────────

HEADER = (
    f"{'configuration':<14}{'cells':>8}{'found':>11}{'allowed':>11}{'no curve yet':>14}{'forbidden':>11}" f"   notes"
)


def result_line(config, summary, complete):
    if summary is None:
        return f"{config:<14}   skipped: no relations yet, since the search for this configuration has not finished"
    notes = []
    if summary["two_row_closures_with_2_to_dim"] is not None:
        notes.append(f"{summary['two_row_closures_with_2_to_dim']}/{summary['cells']} closures have 2^dim cells")
    if not complete:
        notes.append("search unfinished, so 'no curve yet' is too large")
    if summary["found_but_forbidden"]:
        notes.append("CHECK THIS: some found relations are forbidden")
    notes.append(f"{summary['secs']}s")
    return (
        f"{config:<14}{summary['cells']:>8,}{summary['found']:>11,}{summary['allowed']:>11,}"
        f"{summary['predicted_only']:>14,}{summary['found_but_forbidden']:>11,}   {', '.join(notes)}"
    )


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(REPO, "*", "*", "data_*.json*")), key=os.path.getsize)
    wanted = set(sys.argv[1:])
    if wanted:
        files = [f for f in files if config_key(f) in wanted]
        missing = wanted - {config_key(f) for f in files}
        if missing:
            print("No data file found for:", ", ".join(sorted(missing)), flush=True)
    if not files:
        sys.exit("Nothing to do: run this from the repository root, or set REPO=<folder>.")

    print("Rank test. For each configuration, the columns are the number of cells, the relations the search found,")
    print("the pairs the rank test allows, the allowed pairs with no curve yet, and the found relations the test")
    print("forbids. The last column should always be 0.\n")
    print(HEADER, flush=True)

    checked, problems = 0, 0
    context = mp.get_context("spawn" if sys.platform == "darwin" else "fork")
    with context.Pool(WORKERS) as pool:
        for config, summary, complete in pool.imap_unordered(job, files):
            print(result_line(config, summary, complete), flush=True)
            if summary is not None:
                checked += 1
                problems += summary["found_but_forbidden"] > 0

    if problems:
        print(f"\n{problems} of {checked} configurations have found relations that the rank test forbids.")
    else:
        print(
            f"\nChecked {checked} configuration{'s' if checked != 1 else ''}. No found relation breaks the rank test."
        )
