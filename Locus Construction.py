import sympy as sp
from itertools import permutations, combinations, product as iproduct
import numpy as np

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Locus Construction.py"

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
    """Return number ranges for each Jordan block, e.g. '22' -> [[1,2],[3,4]]."""
    digits = [int(d) for d in num_str]
    blocks, start = [], 1
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
    """All permutations that preserve block ordering."""
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
    """Each filling becomes a 0/1 Schubert cell matrix."""
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
    """For each column i, ensure χC[:,i] is in the span of C[:,0..i]."""
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
    """All r-minors for r = 1 .. k-1, using row subset S and first r columns."""
    rows, all_levels = list(range(n)), []
    for r in range(1, k):
        dets = []
        for rs in combinations(rows, r):
            submat = M.extract(list(rs), list(range(r)))
            d = submat[0, 0] if r == 1 else sp.simplify(submat.det())
            dets.append((rs, d))
        all_levels.append((r, dets))
    return all_levels


# ── Closure search (same algorithm as Closure Representation) ─────────────────

def ratio_limit_check(sv, tv, lam, limit_var):
    """lim_{limit_var→∞} (λ·sv[i]) / tv[i] = same nonzero constant for all nonzero tv[i]."""
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
    for (_, src_dets), (_, tgt_dets) in zip(src_levels, tgt_levels):
        sv = [sp.simplify(sp.sympify(d).subs(src_subs)) for _, d in src_dets]
        tv = [sp.simplify(sp.sympify(d).subs(tgt_subs)) for _, d in tgt_dets]
        lam = next(
            (sp.simplify(t / s) for s, t in zip(sv, tv) if t != 0 and s != 0),
            sp.Integer(1)
        )
        if not ratio_limit_check(sv, tv, lam, limit_var):
            return False
    return True


def clean_subs(subs, levels):
    """Remove no-op substitutions: keys not in levels, or values equal to keys."""
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
    """Search for limit_var, src_subs, tgt_subs proving tgt is in src's closure."""
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


# ── Locus / parameter space description ────────────────────────────────────────

def format_plucker_map(levels):
    """Format the Plucker map as a tuple of bracketed vectors, one per level."""
    vecs = ['[' + ' : '.join(str(d) for _, d in dets) + ']' for _, dets in levels]
    if len(vecs) == 1:
        return '{' + vecs[0] + '}'
    return '{' + vecs[0] + ',\n            ' + ',\n            '.join(vecs[1:-1] + [vecs[-1] + '}'])


def format_boundary(boundary_perm, limit_var, src_subs, tgt_subs):
    """One-line description of a boundary stratum in the parameter space."""
    ops = []
    if src_subs:
        ops.append("set " + ", ".join(f"{k}={v}" for k, v in src_subs.items()))
    if limit_var:
        ops.append(f"{limit_var}→∞")
    if tgt_subs:
        ops.append("tgt params: " + ", ".join(f"{k}={v}" for k, v in tgt_subs.items()))
    return f"C_{boundary_perm}:  {';  '.join(ops) if ops else '(direct)'}"


def describe_locus(src_subs, limit_var, free_vars):
    """
    Describe the boundary locus inside A^d for one closure relation.
    Returns a string like "{a4 = 0};  a2 → ∞" or "{a2 = a4**3};  a4 → ∞".
    """
    parts = []
    if src_subs:
        eq_parts = ", ".join(f"{k} = {v}" for k, v in src_subs.items())
        parts.append("{" + eq_parts + "}")

    if limit_var is not None:
        parts.append(f"{limit_var} → ∞")

    # If no limit and no subs, this is just the open interior
    if not parts:
        return "(direct inclusion, no degeneration)"

    return ";  ".join(parts)


def locus_dimension(src_subs, limit_var, free_vars):
    """Estimate the dimension of the boundary locus inside A^d."""
    constrained = set(src_subs.keys()) if src_subs else set()
    if limit_var is not None:
        constrained.add(limit_var)
    free_remaining = [v for v in free_vars if v not in constrained]
    return len(free_remaining)


# ── Main ───────────────────────────────────────────────────────────────────────

num_str = input("Enter number string: ")

digits = [int(d) for d in num_str]
start, parts = 1, []
for d in digits:
    parts.append('(' + ''.join(str(i) for i in range(start, start + d)) + ')')
    start += d
print("Partition: " + ''.join(parts))

tableaux = young_tableaux_from_digits(num_str)
J        = jordan_from_digits(num_str)
valid    = generate_valid_fillings(tableaux, num_str)
cells    = fillings_to_schubert_cells(valid)

# Build Plücker data for every cell
cell_data = []
for idx, cell in enumerate(cells):
    M, a_dict = fill_symbolic_entries(cell)
    C_sym  = sp.Matrix([[a_dict.get(val, val) for val in row] for row in M])
    XC_sym = J * C_sym
    for r in range(J.rows):
        if all(J[r, c] == 0 for c in range(J.cols)):
            for c in range(XC_sym.cols):
                XC_sym[r, c] = 0

    _, free_vars, relations = springer_span_checks(C_sym, XC_sym)
    C_final = C_sym.subs(relations)
    n, k    = C_final.shape
    levels  = plucker_coordinates(C_final, n, k)
    perm    = ''.join(str(v) for row in valid[idx] for v in row)

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
    print(f"  Processed C_{perm}")

# Precompute all closure results
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

# ── Locus descriptions ─────────────────────────────────────────────────────────

print("\n" + "=" * 80)
print("PARAMETER SPACE LOCUS DESCRIPTIONS")
print("=" * 80)

for src in cell_data:
    dim = len(src['tgt_free_vars'])
    if src['tgt_free_vars']:
        coord_str = ", ".join(str(v) for v in src['tgt_free_vars'])
        param_line = f"A^{dim}  with coordinates ({coord_str})"
        map_lhs    = f"({coord_str})"
    else:
        param_line = "A^0  (a single point, no free parameters)"
        map_lhs    = "(*)"

    print(f"\nC_{src['perm']}")
    print(f"  Parameter space  : {param_line}")
    print(f"  Plücker map      : {map_lhs} ↦ {format_plucker_map(src['levels'])}")

    # Find which cells appear as boundary strata in this cell's closure
    boundaries = []
    for tgt in cell_data:
        if src['perm'] == tgt['perm']:
            continue
        matched, lv, ss, ts = all_closure_results[(src['perm'], tgt['perm'])]
        if matched:
            boundaries.append((tgt['perm'], lv, ss, ts))

    if boundaries:
        print(f"  Closure          : A^{dim} ∪ (boundary strata listed below)")
        print("  Boundary loci inside the parameter space:")
        for bperm, lv, ss, ts in boundaries:
            locus = describe_locus(ss, lv, src['tgt_free_vars'])
            ldim  = locus_dimension(ss, lv, src['tgt_free_vars'])
            tgt_info = ""
            if ts:
                tgt_info = "   [tgt params: " + ", ".join(f"{k}={v}" for k, v in ts.items()) + "]"
            print(f"    {locus}   (dim {ldim})  →  C_{bperm}{tgt_info}")
    else:
        print("  Closure          : the cell itself (no boundary; closed orbit)")