"""
Cell closures for a Springer fiber
==================================

What this does
--------------
Consider a nilpotent matrix X whose Jordan blocks have the row lengths of a configuration, such as (3,3,3).
The Springer fiber of X is the set of complete flags V_1 ⊂ ... ⊂ V_n with X·V_j ⊆ V_j for every j. It breaks
into pieces called cells, one for each permutation that keeps every row of the tableau in increasing order.

This script answers one question for every pair of cells: does the closure of the first cell meet the second?
Here, it does so in three steps.

    Step 1. Build every cell. Impose the flag condition on its Schubert cell matrix, then compute its
            Plücker coordinates.
    Step 2. Index every point of every cell where the free variables are 0, 1, -1, T or -T.
    Step 3. For each cell, send a curve toward the boundary and record every cell it lands in.

How to run
----------
    python3 "Computational Geometry.py" 333       configuration (3,3,3)
    python3 "Computational Geometry.py" 1-10      use dashes once a row has 10 or more boxes

The result is saved to ~/data_<configuration>.json. Progress is saved every 30 seconds, so running the same
command again resumes an interrupted run.

Settings (environment variables, all optional)
----------------------------------------------
    WORKERS=8          number of processes (default: every core)
    PRUNE_BRUHAT=0     also search pairs that Bruhat order already rules out
    EXACT_STRINGS=1    print Plücker coordinates exactly as older versions did (slow once n >= 8)

For reviewers
-------------
Part 1 uses SymPy and is a direct translation of the math. Part 2 is the fast search. It replaces the old
pair-by-pair loop with vectorized numpy over thousands of substitutions at once, and it finds exactly the same
first curve for every pair. Part 3 handles checkpoints and worker processes.
"""

import os
import sys
import json
import time
import shutil
import tempfile
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from math import gcd

# Each process does its own linear algebra, and we already run one process per core. As such, every process
# should use a single math thread, or the threads end up competing with each other.
for thread_setting in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(thread_setting, "1")

import numpy as np
import sympy as sp
from flint import fmpq

# ── Settings ─────────────────────────────────────────────────────────────────────────────────────────────────────

CONFIG = sys.argv[1] if len(sys.argv) > 1 else "44"

# Stored in the output file. A checkpoint saved under the old block-by-block labeling is never resumed.
LABELING = "column-first"


def setting(name, default):
    return os.environ.get(name, default)


WORKERS = int(setting("WORKERS", os.cpu_count() or 1))
PRUNE_BRUHAT = setting("PRUNE_BRUHAT", "1") != "0"
EXACT_STRINGS = setting("EXACT_STRINGS", "0") == "1"
CHECKPOINT_SECS = float(setting("CHECKPOINT_SECS", "30"))  # save progress at most this often
CHUNK = int(setting("CHUNK", "8192"))  # substitutions handled in one numpy batch
TERM_BUDGET = int(setting("TERM_BUDGET", "1500000"))  # caps (batch size) x (number of terms) to limit memory
SIGN_LAZY = setting("SIGN_LAZY", "1") != "0"  # only try both signs of a variable once the sign matters
PROGRESS = setting("PROGRESS", "1") != "0"  # print a line per limit variable for large cells


# ══ Part 1. Building the cells ════════════════════════════════════════════════════════════════════════════════════


def row_lengths(config):
    """The row lengths of a configuration. For example, '332' -> [3, 3, 2] and '1-10' -> [1, 10]."""
    parts = config.split("-") if "-" in config else config
    return [int(x) for x in parts]


def labeled_tableau(config):
    """
    The boxes of the tableau, labeled 1, ..., n down each column, starting from the leftmost column.

        33 -> [1][3][5]        321 -> [1][4][6]
              [2][4][6]               [2][5]
                                      [3]
    """
    rows = [[0] * length for length in row_lengths(config)]
    label = 1
    for column in range(max(len(row) for row in rows)):
        for row in rows:
            if column < len(row):
                row[column] = label
                label += 1
    return rows


def jordan_from(config):
    """The nilpotent matrix X. Each box sends its basis vector to the box on its left, and the first column to 0."""
    tableau = labeled_tableau(config)
    n = sum(len(row) for row in tableau)
    X = np.zeros((n, n), dtype=int)
    for row in tableau:
        for column in range(1, len(row)):
            X[row[column - 1] - 1, row[column] - 1] = 1
    return sp.Matrix(X)


def blocks_of(config):
    """The groupings, which are the rows of the tableau. For example, 321 -> [[1, 4, 6], [2, 5], [3]]."""
    return [row[:] for row in labeled_tableau(config)]


def valid_fillings(config):
    """
    Every permutation in which each grouping appears in increasing order, in lexicographic order.
    Each one is returned split into the rows of the tableau.

    In other words, these are the shuffles of the groupings. Taking the smallest available label first gives
    lexicographic order directly, with no need to filter all n! permutations.
    """
    blocks = blocks_of(config)
    lengths = row_lengths(config)
    n = sum(lengths)
    next_in_block = [0] * len(blocks)  # how many labels of each grouping are already used
    sequence, out = [], []

    def extend():
        if len(sequence) == n:
            filling, start = [], 0
            for length in lengths:
                filling.append(list(sequence[start : start + length]))
                start += length
            out.append(filling)
            return
        open_blocks = [b for b in range(len(blocks)) if next_in_block[b] < len(blocks[b])]
        for b in sorted(open_blocks, key=lambda b: blocks[b][next_in_block[b]]):
            sequence.append(blocks[b][next_in_block[b]])
            next_in_block[b] += 1
            extend()
            next_in_block[b] -= 1
            sequence.pop()

    extend()
    return out


def filling_to_cell(filling):
    """The permutation matrix of a filling. Column k has its 1 in row w(k)."""
    word = [label for row in filling for label in row]
    n = len(word)
    matrix = [[0] * n for _ in range(n)]
    for column, label in enumerate(word):
        matrix[label - 1][column] = 1
    return matrix


def fill_symbolic(cell):
    """
    The Schubert cell matrix before the flag condition. Every position with no 1 to its left and no 1 above it
    gets a variable a1, a2, ..., numbered in reading order.
    """
    symbols, index = {}, 1
    matrix = [row[:] for row in cell]
    for r in range(len(matrix)):
        for c in range(len(matrix)):
            if matrix[r][c] != 0:
                continue
            one_to_the_left = any(matrix[r][cc] == 1 for cc in range(c))
            one_above = any(matrix[rr][c] == 1 for rr in range(r))
            if not one_to_the_left and not one_above:
                name = f"a{index}"
                matrix[r][c] = name
                symbols[name] = sp.symbols(name)
                index += 1
    return matrix, symbols


def variable_number(symbol):
    """'a12' -> 12, so variables sort by number rather than alphabetically."""
    return int(symbol.name[1:])


def span_minors(C, XC, j):
    """
    X·v_j lies in the span of v_0, ..., v_j exactly when every (j+2) x (j+2) minor of the matrix
    [v_0 ... v_j | X·v_j] vanishes. Here, we return the minors that are not already zero.

    The minors are built one column at a time by expanding along the last column, so each size reuses the
    minors of the size before it. Rows that are entirely zero are skipped, since any minor using one is 0.
    """
    A = sp.Matrix.hstack(C[:, : j + 1], XC[:, j])
    size = j + 2
    nonzero_rows = [r for r in range(A.rows) if any(A[r, c] != 0 for c in range(size))]
    if len(nonzero_rows) < size:
        return []
    previous = {(r,): sp.expand(A[r, 0]) for r in nonzero_rows}
    for t in range(2, size + 1):
        current = {}
        for rows in combinations(nonzero_rows, t):
            minor = sp.Integer(0)
            for i, row in enumerate(rows):
                entry = A[row, t - 1]
                if entry == 0:
                    continue
                smaller = previous[rows[:i] + rows[i + 1 :]]
                if smaller != 0:
                    minor += (-1) ** (i + t - 1) * entry * smaller
            current[rows] = sp.expand(minor)
        previous = current
    return [previous[rows] for rows in combinations(nonzero_rows, size) if previous[rows] != 0]


def in_span(C, XC, j, pivots):
    """
    A fast exact check that X·v_j lies in the span of v_0, ..., v_j.

    Restricted to their pivot rows, the first j+1 columns form a unit triangular matrix. As such, the
    coefficients come out by forward substitution with no division, and X·v_j is in the span exactly when
    the remainder is zero.
    """
    target = XC[:, j]
    coefficients = []
    for l in range(j + 1):
        known = sum((coefficients[i] * C[pivots[l], i] for i in range(l)), sp.Integer(0))
        coefficients.append(sp.expand(target[pivots[l]] - known))
    for r in range(C.rows):
        fit = sum((coefficients[i] * C[r, i] for i in range(j + 1)), sp.Integer(0))
        if sp.expand(target[r] - fit) != 0:
            return False
    return True


def span_check(C, XC):
    """
    Impose the flag condition X·V_j ⊆ V_j for every j, and return (free variables, relations).

    Each relation is solved for its lowest-numbered variable that appears linearly with a constant coefficient,
    so the remaining free variables are always the same ones.
    """
    all_vars = sorted({s for s in C.free_symbols if s.name.startswith("a")}, key=variable_number)
    relations = {}
    pivots = [
        next(r for r in range(C.rows) if C[r, c] == 1 and all(C[r, cc] == 0 for cc in range(c + 1, C.cols)))
        for c in range(C.cols)
    ]
    for j in range(C.cols):
        while True:
            C_now, XC_now = C.subs(relations), XC.subs(relations)
            if in_span(C_now, XC_now, j, pivots):
                break
            equations = span_minors(C_now, XC_now, j)
            if not equations:
                break
            update = solve_one_linear(equations, all_vars)
            if update is None:
                solutions = sp.solve(equations, [v for v in all_vars if v not in relations], dict=True)
                if not solutions:
                    raise RuntimeError(f"the span condition cannot be satisfied at column {j}")
                print(f"    warning: nonlinear span condition at column {j}, taking branch {solutions[0]}", flush=True)
                update = solutions[0]
            relations = {var: sp.expand(value.subs(update)) for var, value in relations.items()}
            relations.update(update)
    return [v for v in all_vars if v not in relations], relations


def solve_one_linear(equations, all_vars):
    """Solve the simplest equation for a variable that appears linearly, as {variable: value}, or None."""
    for equation in sorted(equations, key=lambda e: (len(sp.Add.make_args(e)), str(e))):
        for var in sorted(equation.free_symbols & set(all_vars), key=variable_number):
            poly = sp.Poly(equation, var)
            if poly.degree() == 1 and not poly.coeffs()[0].free_symbols:
                return {var: sp.expand(sp.solve(equation, var)[0])}
    return None


SIMPLIFIED = {}  # the same minor appears in many cells, so each one is simplified once


def simplified(expression):
    result = SIMPLIFIED.get(expression)
    if result is None:
        result = sp.simplify(expression)
        if len(SIMPLIFIED) > 200000:
            SIMPLIFIED.clear()
        SIMPLIFIED[expression] = result
    return result


def plucker(M, n, k):
    """
    The Plücker coordinates of the flag: for r = 1, ..., k-1, every r x r minor using the first r columns.

    Each level is expanded along its last column and reuses the level before it. Minors that are 0 or a single
    term are already in final form, so only the rest are simplified.
    """
    out, previous = [], {}
    for r in range(1, k):
        minors, current = [], {}
        for rows in combinations(range(n), r):
            if r == 1:
                minor = sp.expand(M[rows[0], 0])
            else:
                minor = sp.Integer(0)
                for t, row in enumerate(rows):
                    entry = M[row, r - 1]
                    if entry == 0:
                        continue
                    smaller = previous[rows[:t] + rows[t + 1 :]]
                    if smaller != 0:
                        minor += (-1) ** (t + r - 1) * entry * smaller
                minor = sp.expand(minor)
            current[rows] = minor
            if r == 1 or not isinstance(minor, sp.Add):
                printed = minor
            elif EXACT_STRINGS:
                printed = sp.simplify(M.extract(list(rows), list(range(r))).det())
            else:
                printed = simplified(minor)
            minors.append((rows, printed))
        out.append((r, minors))
        previous = current
    return out


def compile_levels(levels):
    """
    Rewrite every Plücker coordinate as a plain list of terms, so the search never has to touch SymPy.

        compiled[level][minor] = [(coefficient, {variable: exponent}), ...]
    """
    out = []
    for _, minors in levels:
        entries = []
        for _, minor in minors:
            expression = sp.expand(sp.sympify(minor))
            symbols = sorted(expression.free_symbols, key=variable_number)
            if not symbols:
                terms = [(to_fraction(expression), {})]
            else:
                terms = []
                for monomial, coefficient in sp.Poly(expression, *symbols).terms():
                    exponents = {symbols[i].name: int(e) for i, e in enumerate(monomial) if e}
                    terms.append((to_fraction(coefficient), exponents))
            entries.append(terms)
        out.append(entries)
    return out


def to_fraction(number):
    rational = sp.Rational(number)
    return fmpq(int(rational.p), int(rational.q))


def freeze(compiled):
    """Turn the fractions into plain tuples, so they can be sent to worker processes."""
    return [
        [[((int(c.p), int(c.q)), tuple(sorted(exponents.items()))) for c, exponents in entry] for entry in level]
        for level in compiled
    ]


def build_one_cell(task):
    """The flag condition and Plücker coordinates of one cell. This runs inside a worker process."""
    filling, X = task
    word = [label for row in filling for label in row]
    matrix, symbols = fill_symbolic(filling_to_cell(filling))
    C = sp.Matrix([[symbols.get(entry, entry) for entry in row] for row in matrix])
    free_vars, relations = span_check(C, X * C)
    C_final = C.subs(relations)
    n, k = C_final.shape
    levels = plucker(C_final, n, k)
    used = set()
    for _, minors in levels:
        for _, minor in minors:
            used |= sp.sympify(minor).free_symbols
    separator = "-" if len(word) >= 10 else ""  # labels such as 10 need dashes once n >= 10
    return {
        "perm": separator.join(str(label) for label in word),
        "word": word,
        "levels": levels,
        "frozen": freeze(compile_levels(levels)),
        "free_vars": free_vars,
        "plucker_vars": [v for v in free_vars if v in used],
        "plucker_strs": ["[" + " : ".join(str(m) for _, m in minors) + "]" for _, minors in levels],
    }


def build_cells(config):
    X = jordan_from(config)
    print("Labeled tableau:")
    for row in labeled_tableau(config):
        print("  " + "".join(f"[{label}]" for label in row))
    print("Groupings: " + "".join("(" + "".join(map(str, block)) + ")" for block in blocks_of(config)))
    print("Nilpotent operator X:")
    sp.pprint(X)

    tasks = [(filling, X) for filling in valid_fillings(config)]
    print(f"\nStep 1 of 3: building {len(tasks)} cells", flush=True)
    pool = None
    if WORKERS > 1 and len(tasks) >= 40:
        pool = mp.get_context("spawn").Pool(WORKERS)
        results = pool.imap(build_one_cell, tasks, chunksize=max(1, len(tasks) // (WORKERS * 8)))
    else:
        results = map(build_one_cell, tasks)
    cells = []
    for i, cell in enumerate(results):
        cells.append(cell)
        print(f"  cell {i + 1}/{len(tasks)}: C_{cell['perm']}", flush=True)
    if pool:
        pool.close()
        pool.join()
    return cells


# ══ Part 2. The closure search ════════════════════════════════════════════════════════════════════════════════════
#
# The idea. Take a source cell and substitute 0, ±1 or ±T^k for its free variables, with one variable (the limit
# variable) equal to T. As T goes to infinity, the flag tends to a limit point. If that limit point lies in a
# target cell, then the closure of the source meets the target.
#
# In practice, both points are compared through their Plücker coordinates, level by level and up to scale. Every
# term of a coordinate becomes 0 or ±c·T^d, so the limit of each level is the vector of coefficients of its
# highest power of T. numpy computes these leading vectors for thousands of substitutions at once. Vectors are
# first compared by a 64-bit hash, and every match is then confirmed exactly.
#
# A substitution is stored as three integer arrays, one column per variable:
#     POW   the power of T (0 for the constants)
#     NEG   1 if the value is negative
#     ZERO  1 if the value is 0


def bruhat_le(v, w):
    """v <= w in Bruhat order, where the identity is the smallest."""
    for i in range(1, len(v) + 1):
        if any(a > b for a, b in zip(sorted(v[:i]), sorted(w[:i]))):
            return False
    return True


def clean(substitution, levels):
    """Drop the parts of a substitution that have no effect, before recording it."""
    if not substitution:
        return substitution
    used = set()
    for _, minors in levels:
        for _, minor in minors:
            used |= sp.sympify(minor).free_symbols
    return {var: value for var, value in substitution.items() if var in used and value != var}


# The values a source variable can take, by index: 0, 1, -1, T, T^2, T^3, -T, -T^2, -T^3
MENU_POW = np.array([0, 0, 0, 1, 2, 3, 1, 2, 3], dtype=np.int64)
MENU_NEG = np.array([0, 0, 1, 0, 0, 0, 1, 1, 1], dtype=np.int64)
MENU_ZERO = np.array([1, 0, 0, 0, 0, 0, 0, 0, 0], dtype=np.int64)
MENU_FLIP = np.array([0, 2, 1, 6, 7, 8, 3, 4, 5], dtype=np.int64)  # the index of minus the value
# The values a target variable can take: 0, 1, -1, T, -T (as indices into the menu above)
TARGET_MENU = np.array([0, 1, 2, 3, 6], dtype=np.int64)
# Random odd weights for the 64-bit hash of a vector of leading coefficients
HASH_WEIGHTS = np.random.default_rng(20261015).integers(1, 2**63, size=8192, dtype=np.uint64) | np.uint64(1)


def menu_strings(limit_var):
    """The printed source and target menus, used when recording a curve."""
    if limit_var is None:
        return ["0", "1", "-1"], ["0", "1", "-1"]
    T = sp.Symbol(limit_var)
    source = [0, 1, -1, T, T**2, T**3, -T, -(T**2), -(T**3)]
    target = [0, 1, -1, T, -T]
    return [str(sp.sympify(x)) for x in source], [str(sp.sympify(x)) for x in target]


def digits_of(index, num_vars, base=5):
    """Substitutions are numbered in product order. This turns a number back into one menu choice per variable."""
    if num_vars == 0:
        return np.zeros((len(index), 0), np.int64)
    place = np.array([base ** (num_vars - 1 - m) for m in range(num_vars)], dtype=np.int64)
    return (np.asarray(index, dtype=np.int64)[:, None] // place[None, :]) % base


def hash_part(L, offset):
    """The hash of some columns of the leading vectors, starting at column `offset` of the full vector."""
    return (np.ascontiguousarray(L).view(np.uint64) * HASH_WEIGHTS[offset : offset + L.shape[1]]).sum(axis=1)


def batch_size(terms, level, chunk):
    """How many substitutions to handle at once at this level, so that memory stays within TERM_BUDGET."""
    num_terms = max(1, terms["levels"][level]["E"].shape[0])
    return int(max(64, min(chunk, TERM_BUDGET // num_terms)))


# ── Per-worker state ─────────────────────────────────────────────────────────────────────────────────────────────

WORKER = {}  # each process keeps its own copy of the cells and of everything computed from them


def init_worker(payload, index_dir=None):
    WORKER.clear()
    WORKER["index_dir"] = index_dir
    WORKER["cells"] = payload
    WORKER["terms"] = {}


def term_arrays(cell_index, which):
    """
    The terms of every Plücker coordinate of a cell as numpy arrays, over the variables in `which`
    ('search_vars' for a source cell, 'point_vars' for a target cell). Cached per worker.

    For each level:
        E       one row per term, one column per variable, holding the exponents
        entry   which coordinate each term belongs to
        coef    the integer coefficient of each term
    The float copies let numpy use fast matrix multiplication. Every value in them is a small exact integer.
    """
    key = (cell_index, which)
    if key in WORKER["terms"]:
        return WORKER["terms"][key]
    cell = WORKER["cells"][cell_index]
    variables = cell[which]
    column_of = {var: i for i, var in enumerate(variables)}
    num_vars = len(variables)
    levels = []
    for level in cell["compiled"]:
        # Scale the whole level by a common denominator. This is the same projective point.
        common = 1
        for entry in level:
            for (p, q), _ in entry:
                common = common * q // gcd(common, q)
        entry_of_term, coefficients, exponents = [], [], []
        for entry_index, entry in enumerate(level):
            for (p, q), powers in entry:
                entry_of_term.append(entry_index)
                coefficients.append(p * (common // q))
                row = [0] * num_vars
                for var, power in powers:
                    row[column_of[var]] = power
                exponents.append(row)
        E = np.array(exponents, dtype=np.int64).reshape(len(entry_of_term), num_vars)
        total_per_entry = [0] * len(level)
        for entry_index, c in zip(entry_of_term, coefficients):
            total_per_entry[entry_index] += abs(c)
        if max(total_per_entry, default=0) >= 2**52:
            raise OverflowError("coefficients too large for exact float arithmetic")

        def as_float_transposed(A):
            return np.ascontiguousarray(A.T, dtype=np.float64)

        levels.append(
            {
                "num_entries": len(level),
                "E": E,
                "entry": np.array(entry_of_term, dtype=np.int64),
                "coef": np.array(coefficients, dtype=np.float64),
                "E_T": as_float_transposed(E),
                "E_odd_T": as_float_transposed(E % 2),
                "E_odd": np.ascontiguousarray(E % 2, dtype=np.float64),
                "E_positive_T": as_float_transposed(E > 0),
                "degree_bound": 3 * int(E.sum(axis=1).max(initial=0)) + 1,
            }
        )
    result = {"levels": levels, "num_vars": num_vars}
    WORKER["terms"][key] = result
    return result


def limit_points(terms, POW, NEG, ZERO, lo=0, hi=None, want_signs=False):
    """
    The limit points of a batch of substitutions (one per row), for levels lo, ..., hi-1.

    Returns
        L          the leading vector of each level, side by side, in a canonical form
                   (primitive integer vector with its first nonzero entry positive)
        nonzero    whether every level has a nonzero leading vector
        top        the highest power of T that appears
        signs      (only if want_signs) which variables' signs can change the result
    """
    count = POW.shape[0]
    rows = np.arange(count)
    parts = []
    nonzero = np.ones(count, dtype=bool)
    top_max = np.zeros(count, dtype=np.int64)
    POW_f, NEG_f, ZERO_f = (A.astype(np.float64) for A in (POW, NEG, ZERO))
    signs = np.zeros(POW.shape, dtype=bool) if want_signs else None

    for level in terms["levels"][lo:hi]:
        num_entries, degree_bound = level["num_entries"], level["degree_bound"]
        if level["E"].shape[0] == 0:
            nonzero[:] = False
            parts.append(np.zeros((count, num_entries), np.int64))
            continue
        degree = (POW_f @ level["E_T"]).astype(np.int64)  # power of T in each term
        negative = (NEG_f @ level["E_odd_T"]).astype(np.int64) & 1  # sign of each term
        dead = (ZERO_f @ level["E_positive_T"]) > 0  # terms that contain a zero variable
        value = np.where(dead, 0.0, np.where(negative == 1, -level["coef"], level["coef"]))

        # Usually only the terms of the highest power survive in the limit, so sum just those.
        live_degree = np.where(dead, -1, degree)
        max_degree = live_degree.max(axis=1)
        on_top = (live_degree == max_degree[:, None]) & (max_degree[:, None] >= 0)
        slot = rows[:, None] * num_entries + level["entry"][None, :]
        L = np.rint(
            np.bincount(slot.ravel(), weights=np.where(on_top, value, 0.0).ravel(), minlength=count * num_entries)
        )
        L = L.astype(np.int64).reshape(count, num_entries)
        top = max_degree.copy()
        has = (L != 0).any(axis=1)

        # Rarely, the top terms cancel. Then those rows are redone power by power.
        redo = np.nonzero(~has & (max_degree >= 0))[0]
        if len(redo):
            R = len(redo)
            redo_rows = np.arange(R)
            slot = (redo_rows[:, None] * num_entries + level["entry"][None, :]) * degree_bound + degree[redo]
            by_power = np.rint(
                np.bincount(slot.ravel(), weights=value[redo].ravel(), minlength=R * num_entries * degree_bound)
            )
            by_power = by_power.astype(np.int64).reshape(R, num_entries, degree_bound)
            power_present = (by_power != 0).any(axis=1)
            top_redo = degree_bound - 1 - np.argmax(power_present[:, ::-1], axis=1)
            L[redo] = by_power[redo_rows, :, top_redo]
            top[redo] = top_redo
            has[redo] = power_present.any(axis=1)

        if want_signs:
            # A sign matters when the variable has an odd power in a term of the leading vector.
            matters = (on_top.astype(np.float64) @ level["E_odd"]) > 0
            if len(redo):
                matters[redo] = ((~dead[redo]).astype(np.float64) @ level["E_odd"]) > 0
            signs |= matters

        nonzero &= has
        top_max = np.maximum(top_max, top)
        divisor = np.gcd.reduce(np.abs(L), axis=1)
        divisor[divisor == 0] = 1
        first = L[rows, np.argmax(L != 0, axis=1)]
        parts.append((L // divisor[:, None]) * np.where(first < 0, -1, 1)[:, None])

    L_all = np.ascontiguousarray(np.hstack(parts))
    return (L_all, nonzero, top_max, signs) if want_signs else (L_all, nonzero, top_max)


# ── Level-by-level enumeration ───────────────────────────────────────────────────────────────────────────────────
#
# Level r of the Plücker coordinates only uses the first r columns, so it only depends on some of the variables.
# As such, variables are assigned in the order of the first level they appear in. After each level, a partial
# substitution survives only if its hash so far matches some target that is still open. Variables that appear in
# no coordinate cannot change the limit point, so they are fixed to 0 instead of enumerated.
#
# A partial substitution is stored as a row of menu indices, with -1 for variables not assigned yet.


def level_offsets(terms):
    offsets, total = [], 0
    for level in terms["levels"]:
        offsets.append(total)
        total += level["num_entries"]
    return offsets, total


def first_levels(terms):
    """For each variable, the first level that uses it, or None if no level does."""
    first = [None] * terms["num_vars"]
    for r, level in enumerate(terms["levels"]):
        if level["E"].shape[0] == 0:
            continue
        used = (level["E"] != 0).any(axis=0)
        for var in range(terms["num_vars"]):
            if first[var] is None and used[var]:
                first[var] = r
    return first


def menu_arrays(choices, num_vars, fixed):
    """Turn rows of menu indices (with -1 for unassigned) into the POW, NEG and ZERO arrays."""
    count = choices.shape[0]
    POW = np.zeros((count, num_vars), np.int64)
    NEG = np.zeros_like(POW)
    ZERO = np.ones_like(POW)
    for column, (power, negative, zero) in fixed.items():
        POW[:, column], NEG[:, column], ZERO[:, column] = power, negative, zero
    free = [c for c in range(num_vars) if c not in fixed]
    if free:
        free_choices = choices[:, free]
        assigned = free_choices >= 0
        index = np.where(assigned, free_choices, 0)
        POW[:, free] = MENU_POW[index]
        NEG[:, free] = MENU_NEG[index]
        ZERO[:, free] = np.where(assigned, MENU_ZERO[index], 1)  # an unassigned variable has no effect yet
    return POW, NEG, ZERO


def expand_rows(choices, new_vars, values, start, stop):
    """Rows start, ..., stop-1 of (current rows) x (every choice of values for the new variables)."""
    base = len(values)
    per_row = base ** len(new_vars)
    flat = np.arange(start, stop, dtype=np.int64)
    parent, combination = flat // per_row, flat % per_row
    expanded = choices[parent].copy()
    for m, column in enumerate(new_vars):
        expanded[:, column] = values[(combination // (base ** (len(new_vars) - 1 - m))) % base]
    return expanded, parent


def enumerate_levels(terms, groups, values, fixed, accept, emit, chunk=None):
    """
    Depth-first search, one level at a time. groups[r] lists the variables first used at level r.
    accept(r, hashes, L, nonzero, top) decides which rows survive level r, and emit(choices, hashes) receives
    the rows that survive every level.
    """
    chunk = chunk or CHUNK
    num_vars, num_levels = terms["num_vars"], len(terms["levels"])
    offsets, _ = level_offsets(terms)
    values = np.asarray(values, dtype=np.int64)

    def search(choices, hashes, r):
        if r == num_levels:
            emit(choices, hashes)
            return
        total = choices.shape[0] * len(values) ** len(groups[r])
        step = batch_size(terms, r, chunk)
        for start in range(0, total, step):
            expanded, parent = expand_rows(choices, groups[r], values, start, min(start + step, total))
            L, nonzero, top = limit_points(terms, *menu_arrays(expanded, num_vars, fixed), r, r + 1)
            new_hashes = hashes[parent] + hash_part(L, offsets[r])
            keep = accept(r, new_hashes, L, nonzero, top)
            if keep.any():
                search(expanded[keep], new_hashes[keep], r + 1)

    search(np.full((1, num_vars), -1, np.int64), np.zeros(1, np.uint64), 0)


def enumerate_levels_signed(terms, groups, values, fixed, accept, emit, chunk=None):
    """
    The same search, but new variables first take only nonnegative values. A variable's sign is tried both ways
    only once it can change a leading vector. Until then, + stands for both signs. This is safe because the two
    give the same limit point, and + comes first in product order.
    """
    chunk = chunk or CHUNK
    num_vars, num_levels = terms["num_vars"], len(terms["levels"])
    offsets, _ = level_offsets(terms)
    values = np.asarray(values, dtype=np.int64)

    def one_level(choices, r):
        return limit_points(terms, *menu_arrays(choices, num_vars, fixed), r, r + 1, want_signs=True)

    def search(choices, undecided, hashes, r):  # undecided marks variables whose sign has not been tried yet
        if r == num_levels:
            emit(choices, hashes)
            return
        total = choices.shape[0] * len(values) ** len(groups[r])
        step = max(32, batch_size(terms, r, chunk) // 4)
        for start in range(0, total, step):
            expanded, parent = expand_rows(choices, groups[r], values, start, min(start + step, total))
            undecided_now = undecided[parent].copy()
            for column in groups[r]:
                undecided_now[:, column] = expanded[:, column] != 0
            L, nonzero, _, signs = one_level(expanded, r)
            to_flip = signs & undecided_now
            while to_flip.any():
                # Try the other sign of every variable whose sign now matters.
                first_new = len(expanded)
                for column in np.nonzero(to_flip.any(axis=0))[0]:
                    rows = to_flip[:, column].copy()
                    if not rows.any():
                        continue
                    undecided_now[rows, column] = False
                    to_flip[rows, column] = False
                    flipped = expanded[rows].copy()
                    flipped[:, column] = MENU_FLIP[flipped[:, column]]
                    first_new = min(first_new, len(expanded))
                    expanded = np.vstack([expanded, flipped])
                    undecided_now = np.vstack([undecided_now, undecided_now[rows]])
                    to_flip = np.vstack([to_flip, to_flip[rows]])
                    parent = np.concatenate([parent, parent[rows]])
                    L = np.vstack([L, L[rows]])
                    nonzero = np.concatenate([nonzero, nonzero[rows]])
                if first_new >= len(expanded):
                    break
                # Recompute the flipped rows. Other signs may now matter too.
                L_new, nonzero_new, _, signs_new = one_level(expanded[first_new:], r)
                L[first_new:] = L_new
                nonzero[first_new:] = nonzero_new
                to_flip = np.zeros_like(undecided_now)
                to_flip[first_new:] = signs_new & undecided_now[first_new:]
            new_hashes = hashes[parent] + hash_part(L, offsets[r])
            keep = accept(r, new_hashes, L, nonzero, None)
            if keep.any():
                search(expanded[keep], undecided_now[keep], new_hashes[keep], r + 1)

    search(
        np.full((1, num_vars), -1, np.int64),
        np.zeros((1, num_vars), bool),
        np.zeros(1, np.uint64),
        0,
    )


# ── The index of target points ───────────────────────────────────────────────────────────────────────────────────


def target_entries(cell_index):
    """
    Every finite point of a target cell with its variables in {0, 1, -1, T, -T}.
    Returns the running hash after each level and the number of each substitution.
    """
    cell = WORKER["cells"][cell_index]
    num_vars = len(cell["point_vars"])
    terms = term_arrays(cell_index, "point_vars")
    num_levels = len(terms["levels"])
    first = first_levels(terms)
    groups = [[] for _ in range(num_levels)]
    fixed = {}
    for var in range(num_vars):
        if first[var] is None:
            fixed[var] = (0, 0, 1)
        else:
            groups[first[var]].append(var)

    found = []
    enumerate_levels(
        terms,
        groups,
        TARGET_MENU,
        fixed,
        accept=lambda r, hashes, L, nonzero, top: nonzero & (top == 0),  # finite and nonzero
        emit=lambda choices, hashes: found.append(choices),
    )
    if not found:
        return np.zeros((0, num_levels), np.uint64), np.zeros(0, np.int64)
    choices = np.concatenate(found)

    # Variables in no level still take every value, so the list is expanded for them.
    unused = [var for var in range(num_vars) if first[var] is None]
    if unused:
        repeats = 5 ** len(unused)
        choices = np.repeat(choices, repeats, axis=0)
        extra = digits_of(np.tile(np.arange(repeats), len(choices) // repeats), len(unused))
        for m, var in enumerate(unused):
            choices[:, var] = TARGET_MENU[extra[:, m]]

    position_in_menu = np.zeros(9, np.int64)
    position_in_menu[TARGET_MENU] = np.arange(5)
    place = np.array([5 ** (num_vars - 1 - m) for m in range(num_vars)], dtype=np.int64)
    numbers = (position_in_menu[choices] * place[None, :]).sum(axis=1) if num_vars else np.zeros(len(choices), np.int64)

    L, _, _ = limit_points(terms, MENU_POW[choices], MENU_NEG[choices], MENU_ZERO[choices])
    offsets, _ = level_offsets(terms)
    running = np.zeros((len(choices), num_levels), np.uint64)
    so_far = np.zeros(len(choices), np.uint64)
    for r, level in enumerate(terms["levels"]):
        so_far = so_far + hash_part(L[:, offsets[r] : offsets[r] + level["num_entries"]], offsets[r])
        running[:, r] = so_far
    return running, numbers


def target_entries_job(cell_index):
    return cell_index, target_entries(cell_index)


def finish_index(hashes, owners, numbers):
    """Put every target point together, sorted by (owner cell, substitution number)."""
    num_levels = max((h.shape[1] for h in hashes), default=0)
    hashes = np.concatenate(hashes) if hashes else np.zeros((0, num_levels), np.uint64)
    owners = np.concatenate(owners) if owners else np.zeros(0, np.int32)
    numbers = np.concatenate(numbers) if numbers else np.zeros(0, np.int64)
    order = np.lexsort((numbers, owners))
    WORKER["hashes"], WORKER["owners"], WORKER["numbers"] = hashes[order], owners[order], numbers[order]
    WORKER["owner_start"] = np.searchsorted(WORKER["owners"], np.arange(len(WORKER["cells"]) + 1))


def load_target_index():
    """Load the shared index from disk (memory-mapped), or build it here if there is none."""
    if "hashes" in WORKER:
        return
    folder = WORKER.get("index_dir")
    if folder and os.path.exists(os.path.join(folder, "hashes.npy")):
        for name in ("hashes", "owners", "numbers"):
            WORKER[name] = np.load(os.path.join(folder, name + ".npy"), mmap_mode="r")
        WORKER["owner_start"] = np.load(os.path.join(folder, "owner_start.npy"))
        return
    results = [target_entries(i) for i in range(len(WORKER["cells"]))]
    finish_index(
        [r[0] for r in results],
        [np.full(len(r[1]), i, np.int32) for i, r in enumerate(results)],
        [r[1].astype(np.int64) for r in results],
    )


def build_shared_index(payload, pool=None):
    """Build the index of target points once, and save it to a temporary folder that every worker reads."""
    init_worker(payload)
    count = len(payload)
    results = [None] * count
    if pool is None:
        for i in range(count):
            results[i] = target_entries(i)
    else:
        largest_first = sorted(range(count), key=lambda i: -len(payload[i]["point_vars"]))
        for i, result in pool.map(target_entries_job, largest_first, chunksize=1):
            results[i] = result
    finish_index(
        [r[0] for r in results],
        [np.full(len(r[1]), i, np.int32) for i, r in enumerate(results)],
        [r[1].astype(np.int64) for r in results],
    )
    folder = tempfile.mkdtemp(prefix="springer_index_")
    for name in ("hashes", "owners", "numbers", "owner_start"):
        np.save(os.path.join(folder, name + ".npy"), WORKER[name])
    size = len(WORKER["owners"])
    WORKER.clear()
    return folder, size


def target_point(cell_index, number):
    """The exact limit vector (as bytes) of one target substitution, cached."""
    cache = WORKER.setdefault("points", {})
    key = (cell_index, number)
    if key not in cache:
        num_vars = len(WORKER["cells"][cell_index]["point_vars"])
        choice = TARGET_MENU[digits_of([number], num_vars)]
        L, _, _ = limit_points(
            term_arrays(cell_index, "point_vars"), MENU_POW[choice], MENU_NEG[choice], MENU_ZERO[choice]
        )
        if len(cache) > 200000:
            cache.clear()
        cache[key] = L[0].tobytes()
    return cache[key]


# ── The search from one source cell ──────────────────────────────────────────────────────────────────────────────


def search_source(job):
    """
    The closure search from one source cell to many target cells at once.

    For each target, this returns the first (limit variable, source substitution, target substitution) in
    product order, which is exactly what the original pair-by-pair search found.
    """
    source, targets = job
    load_target_index()
    cells = WORKER["cells"]
    num_cells = len(cells)
    search_vars = cells[source]["search_vars"]
    terms = term_arrays(source, "search_vars")
    num_levels = len(terms["levels"])
    column_of = {var: i for i, var in enumerate(search_vars)}
    first = first_levels(terms)
    owner_start = WORKER["owner_start"]
    target_dims = np.array([len(c["point_vars"]) for c in cells], np.int64)
    open_targets = np.zeros(num_cells, dtype=bool)  # targets with no curve found yet
    for t in targets:
        if len(cells[t]["point_vars"]) <= len(search_vars):  # a closure only meets smaller cells
            open_targets[t] = True
    found = {}

    def index_of_open_targets():
        """The index entries of the open targets, sorted by full hash, plus each level's hashes for pruning."""
        selected = [np.arange(owner_start[t], owner_start[t + 1]) for t in np.nonzero(open_targets)[0]]
        selected = np.concatenate(selected) if selected else np.zeros(0, np.int64)
        hashes = np.asarray(WORKER["hashes"][selected])
        owners = np.asarray(WORKER["owners"][selected])
        numbers = np.asarray(WORKER["numbers"][selected])
        per_level = [np.unique(hashes[:, r]) for r in range(num_levels)]
        order = np.lexsort((numbers, owners, hashes[:, -1])) if num_levels else np.zeros(0, np.int64)
        full = hashes[order, -1] if num_levels else np.zeros(len(selected), np.uint64)
        return per_level, full, owners[order], numbers[order], int(open_targets.sum()), len(selected)

    per_level, full_hash, owner, number, built_for, num_entries = index_of_open_targets()

    for limit_var in [None] + list(search_vars):
        if not open_targets.any() or num_entries == 0:
            break
        if open_targets.sum() <= 0.75 * built_for:  # rebuild once a quarter of the targets are done
            per_level, full_hash, owner, number, built_for, num_entries = index_of_open_targets()
            if num_entries == 0:
                break
        others = [var for var in search_vars if var != limit_var]
        base = 3 if limit_var is None else 9
        k = len(others)
        source_strings, target_strings = menu_strings(limit_var)
        other_columns = [column_of[var] for var in others]
        place = np.array([base ** (k - 1 - m) for m in range(k)], dtype=np.int64)
        groups = [[] for _ in range(num_levels)]
        fixed = {}
        for var in others:
            column = column_of[var]
            if first[column] is None:
                fixed[column] = (0, 0, 1)
            else:
                groups[first[column]].append(column)
        if limit_var is not None:
            fixed[column_of[limit_var]] = (1, 0, 0)  # the limit variable is T

        # A target may not set the limit variable itself to +T, so record where it sits in each target.
        limit_position = np.full(num_cells, -1, np.int64)
        if limit_var is not None:
            for t in np.nonzero(open_targets)[0]:
                if limit_var in cells[t]["point_vars"]:
                    limit_position[t] = cells[t]["point_vars"].index(limit_var)

        def allowed(owners, numbers):
            """Which target substitutions are allowed for this limit variable."""
            ok = np.ones(len(owners), bool)
            dims = target_dims[owners]
            if limit_var is None:  # only 0 and ±1
                rest = numbers.copy()
                for _ in range(int(dims.max(initial=0))):
                    ok &= (rest % 5) < 3
                    rest //= 5
                return ok
            position = limit_position[owners]
            has_limit = position >= 0
            if has_limit.any():
                digit = (numbers[has_limit] // (5 ** (dims[has_limit] - 1 - position[has_limit]))) % 5
                ok[has_limit] = digit != 3
            return ok

        def accept(r, hashes, L, nonzero, top, per_level=per_level):
            known = per_level[r]
            if len(known) == 0:
                return np.zeros(len(hashes), bool)
            i = np.minimum(np.searchsorted(known, hashes), len(known) - 1)
            return nonzero & (known[i] == hashes)

        best = {}  # target -> (source number, target number, source menu choices)

        def emit(choices, hashes):
            lo = np.searchsorted(full_hash, hashes, "left")
            hi = np.searchsorted(full_hash, hashes, "right")
            hits = np.nonzero(hi > lo)[0]
            if len(hits) == 0:
                return
            menu = np.zeros((len(hits), k), np.int64)
            for m, column in enumerate(other_columns):
                if column not in fixed:
                    menu[:, m] = choices[hits, column]
            source_number = (menu * place[None, :]).sum(axis=1) if k else np.zeros(len(hits), np.int64)
            # Every (survivor, index entry) pair with an equal hash
            count = (hi[hits] - lo[hits]).astype(np.int64)
            survivor = np.repeat(np.arange(len(hits)), count)
            starts = lo[hits].astype(np.int64) - np.concatenate([[0], np.cumsum(count)[:-1]])
            entry = np.repeat(starts, count) + np.arange(count.sum())
            owners = owner[entry]
            keep = open_targets[owners]
            survivor, entry, owners = survivor[keep], entry[keep], owners[keep]
            if len(survivor) == 0:
                return
            keep = allowed(owners, number[entry])
            survivor, entry, owners = survivor[keep], entry[keep], owners[keep]
            if len(survivor) == 0:
                return
            # Keep the smallest (source, target) pair for each target.
            key_source, key_target = source_number[survivor], number[entry]
            order = np.lexsort((key_target, key_source, owners))
            firsts = order[np.r_[True, owners[order][1:] != owners[order][:-1]]]
            for f in firsts:
                t = int(owners[f])
                candidate = (int(key_source[f]), int(key_target[f]))
                current = best.get(t)
                if current is None or candidate < current[:2]:
                    best[t] = candidate + (menu[survivor[f]].copy(),)

        started = time.time()
        if SIGN_LAZY:
            values = [0, 1] if limit_var is None else [0, 1, 3, 4, 5]
            enumerate_levels_signed(terms, groups, values, fixed, accept, emit)
        else:
            enumerate_levels(terms, groups, np.arange(base), fixed, accept, emit)

        if PROGRESS and len(search_vars) >= 8:
            step = ([None] + list(search_vars)).index(limit_var) + 1
            print(
                f"    C_{cells[source]['perm']}: limit {limit_var or 'none'} done "
                f"({step}/{len(search_vars) + 1}, {time.time() - started:.0f}s, "
                f"{len(best)} of {int(open_targets.sum())} remaining targets reached)",
                flush=True,
            )

        # Confirm each target's first pair exactly, since the hash alone could collide.
        for t, (_, target_number, menu) in sorted(best.items()):
            POW = np.zeros((1, len(search_vars)), np.int64)
            NEG = np.zeros_like(POW)
            ZERO = np.zeros_like(POW)
            for m, column in enumerate(other_columns):
                POW[0, column] = MENU_POW[menu[m]]
                NEG[0, column] = MENU_NEG[menu[m]]
                ZERO[0, column] = MENU_ZERO[menu[m]]
            if limit_var is not None:
                POW[0, column_of[limit_var]] = 1
            L, _, _ = limit_points(terms, POW, NEG, ZERO)
            if target_point(t, target_number) != L[0].tobytes():
                result = exact_fallback(source, t, limit_var)  # a 64-bit hash collision, which is very rare
                if result is not None:
                    found[t] = result
                    open_targets[t] = False
                continue
            point_vars = cells[t]["point_vars"]
            target_menu = [int(x) for x in digits_of([target_number], len(point_vars))[0]]
            found[t] = (
                limit_var,
                {others[m]: source_strings[int(menu[m])] for m in range(k)},
                {point_vars[m]: target_strings[target_menu[m]] for m in range(len(point_vars))},
            )
            open_targets[t] = False

    return source, targets, found


def exact_fallback(source, target, limit_var):
    """Check one (source, target, limit variable) by brute force. This only runs after a hash collision."""
    cells = WORKER["cells"]
    search_vars = cells[source]["search_vars"]
    terms = term_arrays(source, "search_vars")
    column_of = {var: i for i, var in enumerate(search_vars)}
    others = [var for var in search_vars if var != limit_var]
    base = 3 if limit_var is None else 9
    k = len(others)
    source_strings, target_strings = menu_strings(limit_var)
    point_vars = cells[target]["point_vars"]
    num_vars = len(point_vars)
    limit_position = point_vars.index(limit_var) if (limit_var is not None and limit_var in point_vars) else -1

    points = {}
    lo, hi = WORKER["owner_start"][target], WORKER["owner_start"][target + 1]
    for number in np.sort(np.asarray(WORKER["numbers"][lo:hi])):
        menu = tuple(int(x) for x in digits_of([int(number)], num_vars)[0])
        if limit_var is None and num_vars and max(menu) >= 3:
            continue
        if limit_position >= 0 and menu[limit_position] == 3:
            continue
        points.setdefault(target_point(target, int(number)), menu)

    place = np.array([base ** (k - 1 - m) for m in range(k)], dtype=np.int64)
    total = base**k
    for start in range(0, total, CHUNK):
        index = np.arange(start, min(start + CHUNK, total), dtype=np.int64)
        digits = (index[:, None] // place[None, :]) % base if k else np.zeros((len(index), 0), np.int64)
        POW = np.zeros((len(index), len(search_vars)), np.int64)
        NEG = np.zeros_like(POW)
        ZERO = np.zeros_like(POW)
        for m, var in enumerate(others):
            column = column_of[var]
            POW[:, column] = MENU_POW[digits[:, m]]
            NEG[:, column] = MENU_NEG[digits[:, m]]
            ZERO[:, column] = MENU_ZERO[digits[:, m]]
        if limit_var is not None:
            POW[:, column_of[limit_var]] = 1
        L, nonzero, _ = limit_points(terms, POW, NEG, ZERO)
        for i in np.nonzero(nonzero)[0]:
            target_menu = points.get(L[i].tobytes())
            if target_menu is not None:
                return (
                    limit_var,
                    {others[m]: source_strings[int(digits[i, m])] for m in range(k)},
                    {point_vars[m]: target_strings[target_menu[m]] for m in range(num_vars)},
                )
    return None


# ══ Part 3. Running and saving ════════════════════════════════════════════════════════════════════════════════════


def save_atomic(path, payload):
    """Write to a temporary file first, so an interrupted save never leaves a broken file behind."""
    temporary = path + ".tmp"
    with open(temporary, "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    os.replace(temporary, path)


def load_checkpoint(path):
    """(relations, settled pairs) from an unfinished run, or None if the file is already complete."""
    if not os.path.exists(path):
        return [], set()
    try:
        prior = json.load(open(path))
    except Exception:
        prior = None
    same_run = prior and prior.get("partition") == CONFIG and prior.get("labeling") == LABELING
    if same_run and prior.get("_complete"):
        print(
            f"\ndata_{CONFIG}.json is already complete ({len(prior.get('edges', []))} relations), so there is "
            f"nothing to do. Delete or move it to compute it again.",
            flush=True,
        )
        return None
    if same_run:
        edges, checked = prior.get("edges", []), set(prior.get("_checked", []))
        print(
            f"\nResuming from data_{CONFIG}.json: {len(checked)} pairs already settled, "
            f"{len(edges)} relations found.",
            flush=True,
        )
        return edges, checked
    print(
        f"\ndata_{CONFIG}.json exists but belongs to another configuration or to the old labeling. Starting fresh.",
        flush=True,
    )
    return [], set()


def make_jobs(cells, checked):
    """One job per source cell, listing the targets still to search. The largest sources go first."""
    jobs, pruned = [], 0
    for source, s in enumerate(cells):
        targets = []
        for target, t in enumerate(cells):
            if target == source:
                continue
            key = f"{s['perm']}>{t['perm']}"
            if key in checked:
                continue
            if PRUNE_BRUHAT and not bruhat_le(t["word"], s["word"]):
                checked.add(key)  # Bruhat order rules this pair out
                pruned += 1
                continue
            targets.append(target)
        if targets:
            jobs.append((source, targets))
    jobs.sort(key=lambda job: -len(cells[job[0]]["free_vars"]))
    return jobs, pruned


def new_pool(payload, index_dir=None):
    # Unlike multiprocessing.Pool, ProcessPoolExecutor raises an error if a worker dies (for example, when it runs
    # out of memory) instead of hanging forever.
    return ProcessPoolExecutor(
        WORKERS, mp_context=mp.get_context("spawn"), initializer=init_worker, initargs=(payload, index_dir)
    )


def main():
    print(f"Computing the cell closures for configuration {CONFIG}", flush=True)
    t0 = time.time()
    cells = build_cells(CONFIG)
    N = len(cells)
    print(f"Built {N} cells in {time.time() - t0:.1f}s", flush=True)

    json_cells = [
        {
            "perm": c["perm"],
            "dim": len(c["plucker_vars"]),
            "free": [str(v) for v in c["plucker_vars"]],
            "plucker": c["plucker_strs"],
        }
        for c in cells
    ]
    position_of = {c["perm"]: i for i, c in enumerate(cells)}

    output_path = os.path.join(os.path.expanduser("~"), f"data_{CONFIG}.json")
    state = load_checkpoint(output_path)
    if state is None:
        return
    edges, checked = state

    def output(final=False):
        edges.sort(key=lambda e: (position_of[e["src"]], position_of[e["tgt"]]))
        out = {"partition": CONFIG, "labeling": LABELING, "cells": json_cells, "edges": edges}
        if final:
            out["_complete"] = True
        else:
            out["_checked"] = sorted(checked)
        return out

    save_atomic(output_path, output())

    total_pairs = N * (N - 1)
    jobs, pruned = make_jobs(cells, checked)
    to_search = sum(len(targets) for _, targets in jobs)

    payload = [
        {
            "perm": c["perm"],
            "search_vars": [v.name for v in c["free_vars"]],
            "point_vars": [v.name for v in c["plucker_vars"]],
            "compiled": c["frozen"],
        }
        for c in cells
    ]

    index_dir = None
    if jobs:
        print("\nStep 2 of 3: indexing the target points", flush=True)
        t1 = time.time()
        if WORKERS > 1:
            with new_pool(payload) as index_pool:
                index_dir, num_points = build_shared_index(payload, index_pool)
        else:
            index_dir, num_points = build_shared_index(payload)
        print(f"  {num_points} target points ({time.time() - t1:.1f}s)", flush=True)

    print(
        f"\nStep 3 of 3: searching closures. Of the {total_pairs} ordered pairs, {len(checked)} are settled "
        f"({pruned} ruled out by Bruhat order) and {to_search} remain, using {WORKERS} worker(s).",
        flush=True,
    )

    recorded = {(e["src"], e["tgt"]) for e in edges}  # relations already in the file, so none is saved twice

    def record(source, targets, found):
        """
        Save the results of one source cell. Every relation is written first, and only then are the pairs marked
        as settled. As such, an interruption can never mark a pair as settled without its relation. At worst, the
        resumed run searches that pair again and skips the relation it already has.
        """
        s = cells[source]
        new_edges = []
        for target in targets:
            if target not in found:
                continue
            t = cells[target]
            limit_var, source_subs, target_subs = found[target]
            source_subs = clean({sp.Symbol(k): sp.sympify(v) for k, v in source_subs.items()}, s["levels"])
            target_subs = clean({sp.Symbol(k): sp.sympify(v) for k, v in target_subs.items()}, t["levels"])
            new_edges.append(
                {
                    "src": s["perm"],
                    "tgt": t["perm"],
                    "limit": limit_var,
                    "srcSubs": {str(k): str(v) for k, v in (source_subs or {}).items()},
                    "tgtSubs": {str(k): str(v) for k, v in (target_subs or {}).items()},
                }
            )
        for edge in new_edges:
            if (edge["src"], edge["tgt"]) not in recorded:
                edges.append(edge)
                recorded.add((edge["src"], edge["tgt"]))
        for target in targets:
            checked.add(f"{s['perm']}>{cells[target]['perm']}")

    def remove_index():
        if index_dir:
            shutil.rmtree(index_dir, ignore_errors=True)

    done, last_save, pool = 0, time.time(), None
    try:
        if WORKERS <= 1:
            init_worker(payload, index_dir)
            results = map(search_source, jobs)
        else:
            pool = new_pool(payload, index_dir)
            results = (f.result() for f in as_completed([pool.submit(search_source, job) for job in jobs]))
        for source, targets, found in results:
            record(source, targets, found)
            done += len(targets)
            if time.time() - last_save > CHECKPOINT_SECS:
                save_atomic(output_path, output())
                last_save = time.time()
            print(
                f"  source C_{cells[source]['perm']} done: {done}/{to_search} pairs, "
                f"{len(edges)} relations ({time.time() - t0:.1f}s)",
                flush=True,
            )
        if pool:
            pool.shutdown()
        save_atomic(output_path, output(final=True))
        remove_index()
        print(
            f"\nFinished. Wrote data_{CONFIG}.json with {N} cells and {len(edges)} relations "
            f"in {time.time() - t0:.1f}s.",
            flush=True,
        )
    except (KeyboardInterrupt, Exception) as error:
        worker_died = type(error).__name__ == "BrokenProcessPool"
        if not isinstance(error, KeyboardInterrupt) and not worker_died:
            raise
        if pool:
            pool.shutdown(wait=False, cancel_futures=True)
        if worker_died:
            print("\nA worker process stopped, usually from running out of memory. Try fewer WORKERS.", flush=True)
        save_atomic(output_path, output())
        remove_index()
        print(
            f"\nInterrupted. Progress is saved in data_{CONFIG}.json ({len(checked)}/{total_pairs} pairs settled, "
            f"{len(edges)} relations). Run the same command again to resume.",
            flush=True,
        )
        sys.exit(0)


if __name__ == "__main__":
    main()
