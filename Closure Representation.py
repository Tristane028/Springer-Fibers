import sympy as sp
from itertools import permutations, combinations
import numpy as np

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Closure Representation.py"

# ── Reuse cell-generation machinery ───────────────────────────────────────────

def young_tableaux_from_digits(num_str):
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

def get_block_ranges(num_str):
    digits = [int(d) for d in num_str]
    blocks, start = [], 1
    for d in digits:
        blocks.append(list(range(start, start + d)))
        start += d
    return blocks

def is_order_preserving(perm, blocks):
    perm = list(perm)
    for block in blocks:
        positions = [perm.index(v) for v in block]
        if positions != sorted(positions):
            return False
    return True

def generate_valid_fillings(tableaux, num_str):
    n = sum(len(r) for r in tableaux)
    blocks = get_block_ranges(num_str)
    valid = []
    for perm in permutations(range(1, n+1)):
        if is_order_preserving(perm, blocks):
            m, idx = [], 0
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
    a_symbols, idx = {}, 1
    M = [row[:] for row in cell]
    for r in range(len(M)):
        for c in range(len(M)):
            if M[r][c] == 0:
                if not any(M[r][cc]==1 for cc in range(c)) and \
                   not any(M[rr][c]==1 for rr in range(r)):
                    name = f"a{idx}"
                    M[r][c] = name
                    a_symbols[name] = sp.symbols(name)
                    idx += 1
    return M, a_symbols

def springer_span_checks(C_sym, XC_sym):
    a_vars = sorted(
        {s for s in C_sym.free_symbols if s.name.startswith("a")},
        key=lambda x: int(x.name[1:])
    )
    relations, results = {}, []
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

def plucker_coordinates(M, n, k):
    rows, all_levels = list(range(n)), []
    for r in range(1, k):
        dets = []
        for rs in combinations(rows, r):
            submat = M.extract(list(rs), list(range(r)))
            d = submat[0, 0] if r == 1 else sp.simplify(submat.det())
            dets.append((rs, d))
        all_levels.append((r, dets))
    return all_levels

# ── Closure-check helpers ──────────────────────────────────────────────────────

from itertools import product as iproduct

def ratio_limit_check(sv, tv, lam, limit_var):
    """
    Check lim_{limit_var->inf} (lam*sv[i]) / tv[i] = same nonzero constant for all i.
    Where tv[i]=0, lam*sv[i] must go to 0. If limit_var is None, check exact equality.
    """
    c = None
    for s, t in zip(sv, tv):
        s, t = sp.simplify(s), sp.simplify(t)
        if t == 0:
            val = sp.limit(lam*s, limit_var, sp.oo) if limit_var else sp.simplify(lam*s)
            if sp.simplify(val) != 0:
                return False
        else:
            try:
                ratio = sp.limit(lam*s/t, limit_var, sp.oo) if limit_var                         else sp.simplify(lam*s/t)
                if ratio in (sp.zoo, sp.oo, -sp.oo, 0):
                    return False
                if c is None:
                    c = ratio
                elif sp.simplify(ratio - c) != 0:
                    return False
            except Exception:
                return False
    return c is not None

def try_one_step(src_levels, tgt_levels, limit_var, src_subs, tgt_subs):
    """Apply subs, find uniform λ per level from first nonzero pair, check ratio_limit_check."""
    for (_, src_dets), (_, tgt_dets) in zip(src_levels, tgt_levels):
        sv = [sp.simplify(sp.sympify(d).subs(src_subs)) for _, d in src_dets]
        tv = [sp.simplify(sp.sympify(d).subs(tgt_subs)) for _, d in tgt_dets]
        lam = None
        for s, t in zip(sv, tv):
            if t != 0 and s != 0:
                lam = sp.simplify(t / s)
                break
        if lam is None:
            lam = sp.Integer(1)
        if not ratio_limit_check(sv, tv, lam, limit_var):
            return False
    return True

def is_in_closure(src_levels, tgt_levels, src_free_vars, tgt_free_vars):
    """
    Search over limit_var, src_subs (other src vars -> polynomial in limit_var),
    and tgt_subs (tgt vars -> 0, ±1, ±limit_var).
    Returns (matched, limit_var, src_subs, tgt_subs).
    """
    for limit_var in [None] + list(src_free_vars):
        other_src = [v for v in src_free_vars if v != limit_var]
        if limit_var is not None:
            src_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1),
                        limit_var, limit_var**2, limit_var**3,
                        -limit_var, -limit_var**2, -limit_var**3]
            tgt_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1),
                        limit_var, -limit_var]
        else:
            src_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]
            tgt_vals = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]

        for src_combo in iproduct(src_vals, repeat=len(other_src)):
            src_subs = dict(zip(other_src, src_combo))
            for tgt_combo in iproduct(tgt_vals, repeat=len(tgt_free_vars)):
                tgt_subs = dict(zip(tgt_free_vars, tgt_combo))
                if try_one_step(src_levels, tgt_levels, limit_var, src_subs, tgt_subs):
                    return True, limit_var, src_subs, tgt_subs

    return False, None, None, None

# ── Main ───────────────────────────────────────────────────────────────────────

num_str  = input("Enter number string: ")
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
    for r in range(J.rows):
        if all(J[r, c] == 0 for c in range(J.cols)):
            for c in range(XC_sym.cols):
                XC_sym[r, c] = 0
    _, free_vars, relations = springer_span_checks(C_sym, XC_sym)
    C_final = C_sym.subs(relations)
    n, k    = C_final.shape
    levels  = plucker_coordinates(C_final, n, k)
    perm    = ''.join(str(v) for row in valid[idx] for v in row)
    cell_data.append({
        'perm':      perm,
        'levels':    levels,
        'free_vars': free_vars,
        'tgt_free_vars': free_vars,
    })
    print(f"  Processed C_{perm}")

# ── Closure table ──────────────────────────────────────────────────────────────

print("\n" + "=" * 80)
print("CLOSURE TABLE")
print("Is row cell in the closure of column cell?")
print("=" * 80)

perms = [d['perm'] for d in cell_data]

# Header
col_w = 10
print(f"{'':>{col_w}}", end='')
for p in perms:
    print(f"  C_{p:<6}", end='')
print()
print("-" * (col_w + len(perms) * col_w))

# For each target cell (row), check against each source cell (col)
for tgt in cell_data:
    print(f"  C_{tgt['perm']:<{col_w-3}}", end='')
    for src in cell_data:
        if src['perm'] == tgt['perm']:
            result = ' ✓(self)'
        else:
            matched, _, _, _ = is_in_closure(src['levels'], tgt['levels'],
                                              src['free_vars'], tgt['free_vars'])
            result = '  ✓     ' if matched else '  ✗     '
        print(f"{result:<{col_w}}", end='')
    print()

# ── Detailed closure results ───────────────────────────────────────────────────

def format_limit_vec(src_vec, lam, var):
    """Show λ*src_vec with var→∞ label, using λ to scale all entries uniformly."""
    scaled = [sp.nsimplify(sp.simplify(lam * e)) for e in src_vec]
    entries = ' : '.join(str(e) for e in scaled)
    if var is not None:
        return f"{var}→∞  [{entries}]"
    else:
        return f"[{entries}]"

print("\n" + "=" * 80)
print("CLOSURE DETAILS")
print("=" * 80)

for tgt in cell_data:
    in_closure_of = []
    for src in cell_data:
        if src['perm'] == tgt['perm']:
            continue
        matched, limit_var, src_subs, tgt_subs = is_in_closure(
            src['levels'], tgt['levels'], src['free_vars'], tgt['free_vars'])
        if matched:
            ops = []
            if src_subs:
                ops.append("src: " + ", ".join(f"{k}={v}" for k,v in src_subs.items()))
            if tgt_subs:
                ops.append("tgt: " + ", ".join(f"{k}={v}" for k,v in tgt_subs.items()))
            if limit_var:
                ops.append(f"{limit_var}→∞")
            # Show the scaled source vectors per level
            limit_strs = []
            for i, (_, src_dets) in enumerate(src['levels']):
                sv = [sp.simplify(sp.sympify(d).subs(src_subs or {})) for _,d in src_dets]
                # find lam
                tgt_dets = tgt['levels'][i][1]
                tv = [sp.simplify(sp.sympify(d).subs(tgt_subs or {})) for _,d in tgt_dets]
                lam = None
                for s,t in zip(sv,tv):
                    if t!=0 and s!=0: lam=sp.simplify(t/s); break
                if lam is None: lam=sp.Integer(1)
                limit_strs.append(format_limit_vec(sv, lam, limit_var))
            in_closure_of.append((src['perm'], ops, limit_strs))

    if in_closure_of:
        print(f"\nC_{tgt['perm']} is in the closure of:")
        for src_perm, ops, limit_strs in in_closure_of:
            print(f"  C_{src_perm}  [{', '.join(ops)}]")
            for ls in limit_strs:
                print(f"    {ls}")
    else:
        print(f"\nC_{tgt['perm']} is not in the closure of any other cell.")