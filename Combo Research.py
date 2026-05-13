from itertools import permutations, combinations
import numpy as np
import sympy as sp
from sympy import symbols, Matrix, simplify
from sympy.parsing.sympy_parser import parse_expr

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Combo Research.py"

# ── Young Tableau & Jordan matrix ─────────────────────────────────────────────

def young_tableaux_from_digits(num_str):
    # First digit = bottom row, subsequent digits = rows above
    return [list(range(1, int(d)+1)) for d in num_str]

def jordan_from_digits(num_str):
    digits = [int(d) for d in num_str]
    n = sum(digits)
    J = np.zeros((n, n), dtype=int)
    idx = 0
    for k in digits:
        for i in range(k-1):
            J[idx+i, idx+i+1] = 1
        idx += k
    return sp.Matrix(J)

# ── Valid fillings ─────────────────────────────────────────────────────────────

def is_valid_matrix(matrix):
    for col in range(len(matrix[0])):
        col_vals = [matrix[r][col] for r in range(len(matrix)) if col < len(matrix[r])]
        if sorted(col_vals) != col_vals or len(col_vals) != len(set(col_vals)):
            return False
    return True

def is_column_increasing(matrix):
    max_cols = max(len(r) for r in matrix)
    for c in range(max_cols):
        for r in range(1, len(matrix)):
            if c < len(matrix[r]) and c < len(matrix[r-1]):
                if matrix[r][c] <= matrix[r-1][c]:
                    return False
    return True

def count_inversions(matrix):
    pos = [(val, r, c) for r, row in enumerate(matrix) for c, val in enumerate(row)]
    inv = 0
    for i in range(len(pos)):
        for j in range(i+1, len(pos)):
            vi, ri, ci = pos[i]
            vj, rj, cj = pos[j]
            if ri == rj and ci < cj and vi > vj: inv += 1
            if rj > ri and cj < ci and vj < vi: inv += 1
    return inv

def get_block_ranges(num_str):
    """Return the list of number ranges for each Jordan block."""
    digits = [int(d) for d in num_str]
    blocks = []
    start = 1
    for d in digits:
        blocks.append(list(range(start, start + d)))
        start += d
    return blocks

def is_order_preserving(perm, blocks):
    """Within each block, numbers must appear in increasing order."""
    perm = list(perm)
    for block in blocks:
        positions = [perm.index(v) for v in block]
        if positions != sorted(positions):
            return False
    return True

def generate_valid_fillings(tableaux, num_str):
    n = sum(len(r) for r in tableaux)
    nums = list(range(1, n+1))
    blocks = get_block_ranges(num_str)
    valid = []
    for perm in permutations(nums):
        if is_order_preserving(perm, blocks):
            m = []
            idx = 0
            for row in tableaux:
                m.append(list(perm[idx:idx+len(row)]))
                idx += len(row)
            valid.append(m)
    return valid

def fillings_to_schubert_cells(valid):
    cells = []
    for mat in valid:
        flat = [v for row in mat for v in row]
        n = len(flat)
        S = [[0]*n for _ in range(n)]
        for col, val in enumerate(flat):
            S[val-1][col] = 1
        cells.append(S)
    return cells

def fill_symbolic_entries(cell):
    a_symbols = {}
    idx = 1
    M = [row[:] for row in cell]
    for r in range(len(M)):
        for c in range(len(M)):
            if M[r][c] == 0:
                has_left  = any(M[r][cc] == 1 for cc in range(c))
                has_above = any(M[rr][c] == 1 for rr in range(r))
                if not has_left and not has_above:
                    name = f"a{idx}"
                    M[r][c] = name
                    a_symbols[name] = sp.symbols(name)
                    idx += 1
    return M, a_symbols

# ── Young tableau display ─────────────────────────────────────────────────────

def display_young_tableau(filling):
    """Display a valid filling as a box Young tableau."""
    for row in reversed(filling):
        print("".join(f'[{v}]' for v in row))

# ── Springer span check ────────────────────────────────────────────────────────

def springer_span_checks(C_sym, XC_sym):
    a_vars = sorted(
        {s for s in C_sym.free_symbols if s.name.startswith("a")},
        key=lambda x: int(x.name[1:])
    )
    relations = {}
    results = []

    for i in range(C_sym.cols):
        C_sub  = C_sym.subs(relations)
        XC_sub = XC_sym.subs(relations)
        C_cols = C_sub[:, :i+1]
        target = XC_sub[:, i]

        alphas   = sp.symbols(f"alpha0:{i+1}")
        unknowns = list(alphas) + a_vars
        eqs      = list(C_cols * sp.Matrix(alphas) - target)
        sol      = sp.solve(eqs, unknowns, dict=True)

        if not sol:
            results.append(("NO", {}))
            continue

        s = sol[0]
        new_rel = {}
        for v in a_vars:
            if v in s:
                expr = s[v]
                if any(a in expr.free_symbols for a in alphas):
                    new_rel[v] = "free"
                else:
                    new_rel[v] = expr

        for k, v in new_rel.items():
            if v != "free":
                relations[k] = v

        results.append(("YES", new_rel))

    constrained = set(relations.keys())
    free_vars   = [v for v in a_vars if v not in constrained]
    return results, free_vars, relations

# ── Plücker coordinates ────────────────────────────────────────────────────────

def plucker_coordinates(M, n, k):
    """
    Compute Plücker coordinates from a sympy Matrix M (n x k).
    For each row subset S of size r, takes rows S and the first r columns.
    Runs for r = 1 .. k-1 (excludes the full k x k determinant).
    """
    rows = list(range(n))
    all_levels = []

    for r in range(1, k):
        row_combos = list(combinations(rows, r))
        dets = []
        for rs in row_combos:
            submat = M.extract(list(rs), list(range(r)))
            d = submat[0, 0] if r == 1 else simplify(submat.det())
            dets.append((rs, d))
        all_levels.append((r, dets))

    return all_levels


def print_plucker(all_levels):
    for r, dets in all_levels:
        label = (
            "1-minors (single columns)" if r == 1
            else "2-minors (pairs of columns)" if r == 2
            else f"{r}-minors"
        )
        print(f"\n{label}:")
        coords = []
        for rs, d in dets:
            idx = ''.join(str(i + 1) for i in rs)
            print(f"  p_{idx} = {d}")
            coords.append(d)
        print(f"  [{' : '.join(str(c) for c in coords)}]")

    print("\nPlücker coordinates:")
    parts = ['[' + ' : '.join(str(d) for _, d in dets) + ']'
             for _, dets in all_levels]
    print('(' + ', '.join(parts) + ')')

# ── Main ───────────────────────────────────────────────────────────────────────

num_str  = input("Enter number string: ")
tableaux = young_tableaux_from_digits(num_str)
J        = jordan_from_digits(num_str)

valid = generate_valid_fillings(tableaux, num_str)
cells = fillings_to_schubert_cells(valid)

print("\n=== Schubert Cells and Springer Span Checks ===")

summary = []  # collect per-cell results for the final table

for idx, cell in enumerate(cells, start=1):
    flat = ''.join(str(v) for row in valid[idx-1] for v in row)
    print(f"\n--- Cell {idx} ({flat}) ---")
    print()
    display_young_tableau(valid[idx - 1])
    print()
    for row in cell:
        print(row)

    M, a_dict = fill_symbolic_entries(cell)

    print("\nWith symbolic entries:")
    for row in M:
        print(row)

    C_sym  = sp.Matrix([[a_dict.get(val, val) for val in row] for row in M])
    XC_sym = J * C_sym

    for r in range(J.rows):
        if all(J[r, c] == 0 for c in range(J.cols)):
            for c in range(XC_sym.cols):
                XC_sym[r, c] = 0

    print("\nChi * C:")
    sp.pprint(XC_sym)

    checks, free_vars, relations = springer_span_checks(C_sym, XC_sym)

    print("\nSpringer Span Checks:")
    print("{:<8} {:<8} {:<30}".format("Check", "Y/N", "Relations"))
    print("-" * 40)
    for i, (yn, new_rel) in enumerate(checks, start=1):
        rel = ", ".join(
            f"{k}=free" if v == "free" else f"{k}={v}"
            for k, v in new_rel.items()
        ) if new_rel else "None"
        print("{:<8} {:<8} {:<30}".format(i, yn, rel))

    print("\nFinal Free Vars:", ", ".join(str(v) for v in free_vars))
    print("Final Relations:", ", ".join(f"{k}={v}" for k, v in relations.items()))

    # Apply final relations to C_sym and compute Plücker coordinates
    C_final = C_sym.subs(relations)
    n, k    = C_final.shape

    print("\nPlücker Coordinates:")
    levels = plucker_coordinates(C_final, n, k)
    print_plucker(levels)

    # Collect for summary table
    summary.append({
        "cell":      idx,
        "perm":      ''.join(str(v) for row in valid[idx-1] for v in row),
        "filling":   valid[idx-1],
        "C":         C_final,
        "XC":        J * C_final,
        "free_vars": ', '.join(str(v) for v in free_vars) if free_vars else '—',
        "relations": [f"{k}={v}" for k, v in relations.items()] if relations else ['—'],
        "plucker":   levels,
    })

    print("\n" + "=" * 50)


# ── Summary table ──────────────────────────────────────────────────────────────

def format_plucker_table(levels):
    """Format Plucker levels as a list of strings, one vector per line."""
    vecs = ['[' + ' : '.join(str(d) for _, d in dets) + ']' for _, dets in levels]
    if len(vecs) == 1:
        return ['{' + vecs[0] + '}']
    lines = ['{' + vecs[0] + ',']
    for v in vecs[1:-1]:
        lines.append(' ' + v + ',')
    lines.append(' ' + vecs[-1] + '}')
    return lines

def format_matrix_rows(M):
    """Format a sympy matrix as a list of row strings like [1 0 0 0]."""
    return ['[' + ' '.join(str(M[r, c]) for c in range(M.cols)) + ']'
            for r in range(M.rows)]

def format_tableau(filling):
    """Format a filling as a list of box-row strings, top row first."""
    return [''.join(f'[{v}]' for v in row) for row in reversed(filling)]

print("\n" + "=" * 160)
print("SUMMARY TABLE")
print("=" * 160)

# Fixed column widths
W_CELL  = 6
W_PERM  = 8
W_TAB   = 12
W_C     = 16
W_XC    = 16
W_FREE  = 12
W_REL   = 18

header = (f"{'Perm':<{W_PERM}} {'Tableau':<{W_TAB}} "
          f"{'C':<{W_C}} {'XC':<{W_XC}} {'Free Vars':<{W_FREE}} "
          f"{'Relations':<{W_REL}} Plücker Coordinates")
print(header)
print("-" * 160)

for row in summary:
    tab_lines = format_tableau(row["filling"])
    c_lines   = format_matrix_rows(row["C"])
    xc_lines  = format_matrix_rows(row["XC"])
    plucker_lines = format_plucker_table(row["plucker"])

    n_lines = max(len(tab_lines), len(c_lines), len(xc_lines),
                  len(row["relations"]), len(plucker_lines))
    tab_lines     += [''] * (n_lines - len(tab_lines))
    c_lines       += [''] * (n_lines - len(c_lines))
    xc_lines      += [''] * (n_lines - len(xc_lines))
    plucker_lines += [''] * (n_lines - len(plucker_lines))

    for i in range(n_lines):
        perm_s  = f"C_{row['perm']}" if i == 0 else ''
        free_s  = row["free_vars"]   if i == 0 else ''
        rel_s   = row["relations"][i] if i < len(row["relations"]) else ''
        pluck_s = plucker_lines[i]

        print(f"{perm_s:<{W_PERM}} {tab_lines[i]:<{W_TAB}} "
              f"{c_lines[i]:<{W_C}} {xc_lines[i]:<{W_XC}} {free_s:<{W_FREE}} "
              f"{rel_s:<{W_REL}} {pluck_s}")
    print()