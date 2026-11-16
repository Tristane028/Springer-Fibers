"""
Flag formula
============

What this does
--------------
Prints the matrix of a Springer cell after imposing the flag condition, straight from its permutation.

Previously, this matrix came from solving the span check, which takes a long time once n is large. Here,
nothing is solved. Each entry is read off the permutation w:

    zero       the entry is 0
    free       the entry is a free variable a(r,c), named by its row r and column c
    copied     the entry equals the entry one box to the right in both labels, plus a correction from any
               blocking 1s (see Section 5.4 of the README)

This works for any configuration, including ones with no data file yet.

How to run
----------
    python3 flag_formula.py 333 362951487      configuration (3,3,3), cell 362951487
    python3 flag_formula.py 1-10 1-3-2-...     use dashes once a row has 10 or more boxes

For reviewers
-------------
post_flag_matrix() takes its arithmetic as arguments (var, add, mul, zero, one). As such, the same code
runs symbolically with SymPy here, and with plain integers modulo a prime in the verification scripts.
Every entry is memoized, so each one is computed once.
"""

import sys
from functools import lru_cache

# ── Reading the input ────────────────────────────────────────────────────────────────────────────────────────────


def rows_of(config):
    """The row lengths of a configuration. For example, '332' -> [3, 3, 2] and '1-10' -> [1, 10]."""
    return [int(x) for x in (config.split("-") if "-" in config else config)]


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


def neighbors(rows):
    """L[p] and R[p] are the boxes to the left and right of box p, or None at the ends of a row."""
    labels = [[0] * length for length in rows]
    k = 1
    for column in range(max(rows)):
        for row in labels:
            if column < len(row):
                row[column] = k
                k += 1
    n = sum(rows)
    L, R = [None] * (n + 1), [None] * (n + 1)
    for row in labels:
        for a, b in zip(row, row[1:]):
            R[a], L[b] = b, a
    return L, R


# ── The formula ──────────────────────────────────────────────────────────────────────────────────────────────────


def post_flag_matrix(rows, w, var, add, mul, zero, one, sub=None):
    """
    The matrix C of the cell w after imposing the flag, as a dict {(row, column): entry}, counting from 1.

    var(r, c) creates the free variable at row r, column c. add, mul, zero and one give the arithmetic.
    """
    n = len(w)
    sub = sub or (lambda a, b: add(a, mul(-1, b)))
    L, R = neighbors(rows)
    w = [None] + list(w)  # w[c] is the row of the 1 in column c
    col = [None] * (n + 1)  # col[r] is the column whose 1 is in row r
    for c in range(1, n + 1):
        col[w[c]] = c

    @lru_cache(maxsize=None)
    def entry(r, c):
        q = w[c]
        if r == q:
            return one
        if r > q or col[r] < c:  # below the 1 of this column, or right of the 1 of this row
            return zero
        if R[q] is None:  # q ends its row, so the whole column is free
            return var(r, c)
        k = col[R[q]]  # the column of the box right of q
        if col[r] < k:  # the 1 of row r sits between the two columns, so a new free variable appears
            return var(r, c)
        return sub(shifted(r, k), correction(r, c, k, col[r]))  # copied, minus the correction from blockers

    def shifted(r, k):
        """(X C_k)[r], which is C[R(r), k]. In other words, column k slid up by one box."""
        return entry(R[r], k) if R[r] is not None else zero

    def correction(r, kp, k, upto):
        """The sum of gamma(kp, c) * C[r, c] over the blocking columns c < upto whose 1 lies between r and q."""
        q = w[kp]
        total = zero
        for c in range(1, upto):
            if c != kp and c < k and r < w[c] < q:
                e = entry(r, c)
                if e != zero:
                    total = add(total, mul(gamma(kp, c), e))
        return total

    @lru_cache(maxsize=None)
    def gamma(kp, c):
        """How much of column c is needed to make X C_k equal C_kp plus earlier columns (the cascade lemma)."""
        q, k, r = w[kp], col[R[w[kp]]], w[c]
        g = sub(shifted(r, k), correction(r, kp, k, c))
        if kp < c:  # the 1 of column c sits between the two columns
            g = sub(g, entry(r, kp))
        return g

    return {(r, c): entry(r, c) for r in range(1, n + 1) for c in range(1, n + 1)}


def symbolic(config, perm):
    """The matrix as a SymPy matrix, with variables named a_r_c."""
    import sympy as sp

    rows, w = rows_of(config), parse_perm(perm)
    C = post_flag_matrix(
        rows,
        w,
        var=lambda r, c: sp.Symbol(f"a_{r}_{c}"),
        add=lambda a, b: a + b,
        mul=lambda a, b: sp.expand(a * b),
        zero=sp.Integer(0),
        one=sp.Integer(1),
    )
    n = len(w)
    return sp.Matrix(n, n, lambda i, j: sp.expand(C[(i + 1, j + 1)]))


# ── Printing ─────────────────────────────────────────────────────────────────────────────────────────────────────


def pretty(expression):
    """Write a_7_4 as a(7,4), powers as ^, and drop the * signs, so the matrix is easier to read."""
    import re

    text = str(expression).replace("**", "^").replace("*", "")
    return re.sub(r"a_(\d+)_(\d+)", r"a(\1,\2)", text)


def position(symbol):
    return tuple(map(int, symbol.name.split("_")[1:]))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python3 flag_formula.py <configuration> <cell>, for example 333 362951487")
    config, perm = sys.argv[1], sys.argv[2]
    M = symbolic(config, perm)
    free = sorted(M.free_symbols, key=position)

    shape = ",".join(map(str, rows_of(config)))
    print(f"C_{perm} in configuration ({shape})")
    print(f"Dimension {len(free)}. Free variables: {', '.join(pretty(v) for v in free)}\n")

    # Aligned table with row and column numbers
    table = [[pretty(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    widths = [max(len(str(j + 1)), max(len(table[i][j]) for i in range(M.rows))) for j in range(M.cols)]
    label_width = len(str(M.rows))
    print(" " * (label_width + 3) + "  ".join(str(j + 1).rjust(widths[j]) for j in range(M.cols)))
    for i in range(M.rows):
        cells = "  ".join(table[i][j].rjust(widths[j]) for j in range(M.cols))
        print(f"{str(i + 1).rjust(label_width)} | {cells}")

    polynomials = [(i + 1, j + 1, M[i, j]) for i in range(M.rows) for j in range(M.cols) if len(M[i, j].args) > 1]
    polynomials = [(r, c, e) for r, c, e in polynomials if e.is_Add]
    if polynomials:
        print("\nEntries that are polynomials (a copy plus a correction from blocking 1s):")
        for r, c, e in polynomials:
            print(f"  row {r}, column {c}: {pretty(e)}")
    else:
        print("\nEvery entry is 0, 1 or a single variable, so nothing is corrected in this cell.")
