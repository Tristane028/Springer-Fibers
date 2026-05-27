import sympy as sp
from itertools import permutations, combinations, product as iproduct
import numpy as np

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Closure Representation.py"

# ── Young tableau & Jordan matrix ─────────────────────────────────────────────

def young_tableaux_from_digits(num_str):
    """Build tableau shape: first digit = bottom row, each next digit = row above."""
    return [list(range(1, int(d) + 1)) for d in num_str]


def jordan_from_digits(num_str):
    """Build the nilpotent Jordan matrix for the given partition."""
    digits = [int(d) for d in num_str]
    n = sum(digits)
    J = np.zeros((n, n), dtype=int)
    idx = 0
    for k in digits:
        for i in range(k - 1):
            J[idx + i, idx + i + 1] = 1
        idx += k
    return sp.Matrix(J)


# ── Valid fillings ─────────────────────────────────────────────────────────────

def get_block_ranges(num_str):
    """Return the number ranges for each Jordan block, e.g. '22' -> [[1,2],[3,4]]."""
    digits = [int(d) for d in num_str]
    blocks, start = [], 1
    for d in digits:
        blocks.append(list(range(start, start + d)))
        start += d
    return blocks


def is_order_preserving(perm, blocks):
    """Check that within each block, numbers appear in increasing order."""
    perm = list(perm)
    for block in blocks:
        positions = [perm.index(v) for v in block]
        if positions != sorted(positions):
            return False
    return True


def generate_valid_fillings(tableaux, num_str):
    """Generate all permutations that preserve block ordering."""
    n = sum(len(r) for r in tableaux)
    blocks = get_block_ranges(num_str)
    valid = []
    for perm in permutations(range(1, n + 1)):
        if not is_order_preserving(perm, blocks):
            continue
        m, idx = [], 0
        for row in tableaux:
            m.append(list(perm[idx:idx + len(row)]))
            idx += len(row)
        valid.append(m)
    return valid


def fillings_to_schubert_cells(valid):
    """Convert each filling to a 0/1 Schubert cell matrix."""
    cells = []
    for mat in valid:
        flat = [v for row in mat for v in row]
        n = len(flat)
        S = [[0] * n for _ in range(n)]
        for col, val in enumerate(flat):
            S[val - 1][col] = 1
        cells.append(S)
    return cells


def fill_symbolic_entries(cell):
    """Insert free symbolic variables a1, a2, ... where no death ray blocks the entry."""
    a_symbols, idx = {}, 1
    M = [row[:] for row in cell]
    for r in range(len(M)):
        for c in range(len(M)):
            if M[r][c] != 0:
                continue
            has_left  = any(M[r][cc] == 1 for cc in range(c))
            has_above = any(M[rr][c] == 1 for rr in range(r))
            if not has_left and not has_above:
                name = f"a{idx}"
                M[r][c] = name
                a_symbols[name] = sp.symbols(name)
                idx += 1
    return M, a_symbols


# ── Springer span check ────────────────────────────────────────────────────────

def springer_span_checks(C_sym, XC_sym):
    """
    For each column i, check that χC[:,i] lies in span(C[:,0..i]).
    Returns (results, free_vars, relations) where relations maps constrained
    a-variables to their values.
    """
    a_vars = sorted(
        {s for s in C_sym.free_symbols if s.name.startswith("a")},
        key=lambda x: int(x.name[1:])
    )
    relations, results = {}, []

    for i in range(C_sym.cols):
        C_sub  = C_sym.subs(relations)
        XC_sub = XC_sym.subs(relations)
        C_cols = C_sub[:, :i + 1]
        target = XC_sub[:, i]

        alphas   = sp.symbols(f"alpha0:{i + 1}")
        unknowns = list(alphas) + a_vars
        eqs      = list(C_cols * sp.Matrix(alphas) - target)
        sol      = sp.solve(eqs, unknowns, dict=True)

        if not sol:
            results.append(("NO", {}))
            continue

        s, new_rel = sol[0], {}
        for v in a_vars:
            if v in s:
                expr = s[v]
                new_rel[v] = "free" if any(a in expr.free_symbols for a in alphas) else expr

        for k, v in new_rel.items():
            if v != "free":
                relations[k] = v

        results.append(("YES", new_rel))

    free_vars = [v for v in a_vars if v not in relations]
    return results, free_vars, relations


# ── Plücker coordinates ────────────────────────────────────────────────────────

def plucker_coordinates(M, n, k):
    """
    Compute all r-minors for r = 1 .. k-1.
    For row subset S of size r: det of rows S × first r columns.
    """
    rows, all_levels = list(range(n)), []
    for r in range(1, k):
        dets = []
        for rs in combinations(rows, r):
            submat = M.extract(list(rs), list(range(r)))
            d = submat[0, 0] if r == 1 else sp.simplify(submat.det())
            dets.append((rs, d))
        all_levels.append((r, dets))
    return all_levels


# ── Closure helpers ────────────────────────────────────────────────────────────

def ratio_limit_check(sv, tv, lam, limit_var):
    """
    Verify that lim_{limit_var → ∞} (λ·sv[i]) / tv[i] equals the same nonzero
    constant for every i where tv[i] ≠ 0, and that λ·sv[i] → 0 where tv[i] = 0.
    If limit_var is None, performs exact equality check instead.
    """
    constant = None
    for s, t in zip(sv, tv):
        s, t = sp.simplify(s), sp.simplify(t)
        if t == 0:
            val = sp.limit(lam * s, limit_var, sp.oo) if limit_var else sp.simplify(lam * s)
            if sp.simplify(val) != 0:
                return False
        else:
            try:
                ratio = (sp.limit(lam * s / t, limit_var, sp.oo)
                         if limit_var else sp.simplify(lam * s / t))
                if ratio in (sp.zoo, sp.oo, -sp.oo, 0):
                    return False
                if constant is None:
                    constant = ratio
                elif sp.simplify(ratio - constant) != 0:
                    return False
            except Exception:
                return False
    return constant is not None


def try_one_step(src_levels, tgt_levels, limit_var, src_subs, tgt_subs):
    """
    Apply substitutions to source and target, find a uniform scalar λ per level
    from the first nonzero (target / source) pair, then verify ratio_limit_check.
    """
    for (_, src_dets), (_, tgt_dets) in zip(src_levels, tgt_levels):
        sv = [sp.simplify(sp.sympify(d).subs(src_subs)) for _, d in src_dets]
        tv = [sp.simplify(sp.sympify(d).subs(tgt_subs)) for _, d in tgt_dets]

        # Determine λ from the first nonzero (target, source) pair
        lam = next(
            (sp.simplify(t / s) for s, t in zip(sv, tv) if t != 0 and s != 0),
            sp.Integer(1)
        )

        if not ratio_limit_check(sv, tv, lam, limit_var):
            return False
    return True


def clean_subs(subs, levels):
    """Remove substitutions that have no effect on the given Plucker levels:
    either the key does not appear in the levels, or the value equals the key (identity)."""
    if not subs:
        return subs
    active_syms = {
        sym
        for _, dets in levels
        for _, d in dets
        for sym in sp.sympify(d).free_symbols
    }
    return {
        k: v for k, v in subs.items()
        if k in active_syms and sp.simplify(v - k) != 0
    }


def is_in_closure(src_levels, tgt_levels, src_free_vars, tgt_free_vars):
    """
    Search for a proof that tgt is in the closure of src by finding:
      - limit_var : one src free variable to send to ∞ (or None)
      - src_subs  : map other src free vars to polynomials in limit_var (degree ≤ 3)
      - tgt_subs  : map tgt free vars to constants or ± limit_var

    Returns (matched, limit_var, src_subs, tgt_subs).
    """
    for limit_var in [None] + list(src_free_vars):
        other_src = [v for v in src_free_vars if v != limit_var]

        if limit_var is not None:
            src_vals = [
                sp.Integer(0), sp.Integer(1), sp.Integer(-1),
                limit_var,  limit_var**2,  limit_var**3,
                -limit_var, -limit_var**2, -limit_var**3,
            ]
            tgt_vals = [
                sp.Integer(0), sp.Integer(1), sp.Integer(-1),
                limit_var, -limit_var,
            ]
        else:
            src_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]
            tgt_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]

        for src_combo in iproduct(src_vals, repeat=len(other_src)):
            src_subs = dict(zip(other_src, src_combo))
            for tgt_combo in iproduct(tgt_vals, repeat=len(tgt_free_vars)):
                tgt_subs = dict(zip(tgt_free_vars, tgt_combo))
                # Skip if any tgt substitution is an identity (v = k)
                if any(sp.simplify(v - k) == 0 for k, v in tgt_subs.items()):
                    continue
                if try_one_step(src_levels, tgt_levels, limit_var, src_subs, tgt_subs):
                    return (
                        True,
                        limit_var,
                        clean_subs(src_subs, src_levels),
                        clean_subs(tgt_subs, tgt_levels),
                    )

    return False, None, None, None


# ── Main ───────────────────────────────────────────────────────────────────────

num_str  = input("Enter number string: ")

# Print partition notation: "32" -> (123)(45), "221" -> (12)(34)(5)
digits = [int(d) for d in num_str]
start, parts = 1, []
for d in digits:
    parts.append('(' + ''.join(str(i) for i in range(start, start + d)) + ')')
    start += d
print("\nPartition: " + ''.join(parts))

tableaux = young_tableaux_from_digits(num_str)
J        = jordan_from_digits(num_str)
valid    = generate_valid_fillings(tableaux, num_str)
cells    = fillings_to_schubert_cells(valid)

# Build Plücker data for each cell
cell_data = []
for idx, cell in enumerate(cells):
    M, a_dict = fill_symbolic_entries(cell)
    C_sym  = sp.Matrix([[a_dict.get(val, val) for val in row] for row in M])
    XC_sym = J * C_sym

    # Zero out rows where J acts as zero
    for r in range(J.rows):
        if all(J[r, c] == 0 for c in range(J.cols)):
            for c in range(XC_sym.cols):
                XC_sym[r, c] = 0

    _, free_vars, relations = springer_span_checks(C_sym, XC_sym)
    C_final = C_sym.subs(relations)
    n, k    = C_final.shape
    levels  = plucker_coordinates(C_final, n, k)
    perm    = ''.join(str(v) for row in valid[idx] for v in row)

    # Only track free vars that actually appear in the Plücker coordinates
    plucker_symbols = {
        sym
        for _, dets in levels
        for _, d in dets
        for sym in sp.sympify(d).free_symbols
    }
    active_free = [v for v in free_vars if v in plucker_symbols]

    cell_data.append({
        'perm':          perm,
        'levels':        levels,
        'free_vars':     free_vars,
        'tgt_free_vars': active_free,
    })
    print(f"\n  Processed C_{perm}")

# Precompute all closure results (reused across table, details, and summary)
print("\nComputing closures...")
all_closure_results = {
    (src['perm'], tgt['perm']): is_in_closure(
        src['levels'], tgt['levels'],
        src['free_vars'], tgt['tgt_free_vars']
    )
    for src in cell_data
    for tgt in cell_data
    if src['perm'] != tgt['perm']
}

# ── Closure table ──────────────────────────────────────────────────────────────

print("\n" + "=" * 80)
print("CLOSURE TABLE  —  ✓ = row cell is in the closure of column cell")
print("=" * 80)

perms = [d['perm'] for d in cell_data]
col_w = 10

print(f"{'':>{col_w}}", end='')
for p in perms:
    print(f"  C_{p:<6}", end='')
print()
print("-" * (col_w + len(perms) * col_w))

for tgt in cell_data:
    print(f"  C_{tgt['perm']:<{col_w - 3}}", end='')
    for src in cell_data:
        if src['perm'] == tgt['perm']:
            symbol = ' ✓(self)'
        else:
            matched = all_closure_results[(src['perm'], tgt['perm'])][0]
            symbol  = '  ✓     ' if matched else '  ✗     '
        print(f"{symbol:<{col_w}}", end='')
    print()

# ── Closure details ────────────────────────────────────────────────────────────

def format_scaled_level(src_dets, src_subs, tgt_dets, tgt_subs, limit_var):
    """
    Show λ·(src after src_subs) for one Plücker level, where λ is determined
    from the first nonzero (target / source) pair.
    """
    sv = [sp.simplify(sp.sympify(d).subs(src_subs or {})) for _, d in src_dets]
    tv = [sp.simplify(sp.sympify(d).subs(tgt_subs or {})) for _, d in tgt_dets]

    lam = next(
        (sp.simplify(t / s) for s, t in zip(sv, tv) if t != 0 and s != 0),
        sp.Integer(1)
    )
    scaled  = [sp.nsimplify(sp.simplify(lam * e)) for e in sv]
    entries = ' : '.join(str(e) for e in scaled)

    return f"{limit_var}→∞  [{entries}]" if limit_var else f"[{entries}]"


print("\n" + "=" * 80)
print("CLOSURE DETAILS")
print("=" * 80)

for tgt in cell_data:
    in_closure_of = []
    for src in cell_data:
        if src['perm'] == tgt['perm']:
            continue
        matched, limit_var, src_subs, tgt_subs = all_closure_results[
            (src['perm'], tgt['perm'])
        ]
        if not matched:
            continue

        # Build operation label
        ops = []
        if src_subs:
            ops.append("src: " + ", ".join(f"{k}={v}" for k, v in src_subs.items()))
        if tgt_subs:
            ops.append("tgt: " + ", ".join(f"{k}={v}" for k, v in tgt_subs.items()))
        if limit_var:
            ops.append(f"{limit_var}→∞")

        level_strs = [
            format_scaled_level(
                src['levels'][i][1], src_subs,
                tgt['levels'][i][1], tgt_subs,
                limit_var
            )
            for i in range(len(src['levels']))
        ]
        in_closure_of.append((src['perm'], ops, level_strs))

    if in_closure_of:
        print(f"\nC_{tgt['perm']} is in the closure of:")
        for src_perm, ops, level_strs in in_closure_of:
            print(f"  C_{src_perm}  [{', '.join(ops)}]")
            for ls in level_strs:
                print(f"    {ls}")
    else:
        print(f"\nC_{tgt['perm']} is not in the closure of any other cell.")

# ── Closure summary ────────────────────────────────────────────────────────────

closure_map = {d['perm']: [] for d in cell_data}
for src in cell_data:
    for tgt in cell_data:
        if src['perm'] != tgt['perm']:
            if all_closure_results[(src['perm'], tgt['perm'])][0]:
                closure_map[src['perm']].append(tgt['perm'])

print("\n" + "=" * 80)
print("CLOSURE SUMMARY  —  cells contained in the closure of each permutation")
print("=" * 80)
print(f"{'Permutation':<16}  Cells in closure")
print("-" * 80)
for src in cell_data:
    p       = src['perm']
    members = closure_map[p]
    member_str = ",  ".join(f"C_{m}" for m in members) if members else "(none)"
    print(f"C_{p:<13}  {member_str}")
