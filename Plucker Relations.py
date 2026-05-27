import sympy as sp
from itertools import permutations, combinations, product as iproduct
import numpy as np

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Plucker Relations.py"

# ── Young tableau & Jordan matrix ─────────────────────────────────────────────

def young_tableaux_from_digits(num_str):
    return [list(range(1, int(d) + 1)) for d in num_str]


def jordan_from_digits(num_str):
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


def springer_span_checks(C_sym, XC_sym):
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


# ── Plücker relations ─────────────────────────────────────────────────────────

def plucker_symbol(level, idx_tuple):
    """Build a sympy symbol like p^2_{12} for the 2-Plücker with indices (1,2)."""
    return sp.Symbol("p" + str(level) + "_" + "".join(str(x) for x in idx_tuple))


def plucker_relations_in_level(n, k):
    r"""
    Quadratic Plücker relations for G(k, n), indexed by (k+2)-subsets of {1,..,n}.

    For each S = {s_1 < ... < s_{k+2}}, split into I = first k-1 elements,
    J = last k+1 elements, and write the standard alternating sum:
        sum_{s=0}^{k} (-1)^s p_{I ∪ {j_s}} p_{J \ {j_s}} = 0.

    For G(2,n) this returns the C(n,4) three-term relations.
    """
    if k <= 1 or k >= n:
        return []
    rels = []
    for S in combinations(range(1, n + 1), k + 2):
        I = S[:k - 1]
        J = S[k - 1:]
        I_set = set(I)
        terms = []
        for s, js in enumerate(J):
            I_aug = tuple(sorted(I_set | {js}))
            J_red = tuple(j for j in J if j != js)
            sign  = (-1) ** s
            terms.append(sign * plucker_symbol(k, I_aug) * plucker_symbol(k, J_red))
        rel = sp.expand(sp.Add(*terms))
        if rel != 0:
            rels.append((rel, f"S={S}"))
    return rels


def incidence_relations(n, r):
    r"""
    Incidence relations between level r and level r+1 (V_r ⊂ V_{r+1}).
    For each (r+2)-subset T = {t_1 < ... < t_{r+2}} of {1,...,n}:
        sum_{s=0}^{r+1} (-1)^s p^r_{T \ {t_s}}_short * p^{r+1}_{T \ {t_s}}_long = 0
    
    Concretely for r=1 (linking 1-Plückers to 2-Plückers):
        for each {i<j<k}:  p^1_i p^2_{jk} - p^1_j p^2_{ik} + p^1_k p^2_{ij} = 0
    """
    rels = []
    for T in combinations(range(1, n + 1), r + 2):
        terms = []
        for s, ts in enumerate(T):
            remaining = tuple(x for x in T if x != ts)
            # 1-Plücker on the single index ts; (r+1)-Plücker on the remaining
            # For general r, p^r is on a r-subset and p^{r+1} on an (r+1)-subset.
            # Using the cofactor expansion: pick ts to go to the r-side as a 1-element,
            # but more generally we pick which element goes to which.
            # The clean rule: for an (r+2)-subset T, expand the determinant of the
            # composite matrix [r+1 rows; r rows] on columns T:
            #   sum_s (-1)^s * p^r_{T\{ts}, except first r entries} * p^{r+1}_{T \ {ts}}
            # This is exactly the same as the Plücker formula but mixing levels.
            # For r=1: p^1 has size 1, p^{r+1}=p^2 has size 2, T has size 3.
            # For each ts in T:  (-1)^s * p^1_{ts} * p^{r+1}_{T\{ts}}
            sign = (-1) ** s
            p_low  = plucker_symbol(r,     (ts,) if r == 1 else None)
            # For general r we'd need a different split; here we handle the common r=1 case
            # explicitly and treat higher r as the general expansion.
            if r == 1:
                p_low  = plucker_symbol(1, (ts,))
                p_high = plucker_symbol(2, remaining)
                terms.append(sign * p_low * p_high)
            else:
                # General: split T into a single element ts (going to level r as r-th)
                # is not quite right. Use the cofactor expansion:
                # det of (r+2)×(r+2) matrix where top r+1 rows are M_{r+1} and bottom is M_r.
                # Expansion along the bottom row:
                #   sum_s (-1)^{r+1+s} M_r[col_ts] * det(M_{r+1} on cols T\{ts})
                # But M_r[col_ts] is just the column vector, not a single number...
                # For r > 1, M_r is r×n so picking one column gives an r-vector, not a scalar.
                # The Plücker p^r on a single index doesn't make sense for r > 1.
                #
                # So incidence relations between r and r+1 for r > 1 are more complex.
                # We'd need: for r-subset A and (r+1)-subset B of T (|T|=r+2):
                #   sum over choices...
                # 
                # For now, implement only r=1 ↔ r=2 incidence relations, which is the
                # main case for 2-row Springer fibers (k=2 means we have levels 1,2,3).
                pass
        if r == 1 and terms:
            rel = sp.expand(sp.Add(*terms))
            if rel != 0:
                rels.append((rel, f"T={T}"))
    return rels


# ── Verification ───────────────────────────────────────────────────────────────

def verify_plucker_relations(cell, n):
    """
    Substitute each cell's computed Plücker coordinates into the universal
    relations and check they evaluate to zero.

    Returns a list of (level_or_'inc', relation_polynomial, evaluated_result) for
    each non-trivially-zero verification (i.e., relations where the substitution
    is not the trivial 0 = 0 case before simplification).
    """
    # Build substitution map from symbolic p^r_I to the actual polynomial in a's
    subs_map = {}
    for r, dets in cell['levels']:
        for idx_tuple, expr in dets:
            sym = plucker_symbol(r, tuple(i + 1 for i in idx_tuple))
            subs_map[sym] = expr

    results = []

    # Intra-level relations
    for r, _ in cell['levels']:
        for rel, label in plucker_relations_in_level(n, r):
            substituted = sp.expand(rel.subs(subs_map))
            # Show the substituted form before simplification
            results.append(('level', r, rel, substituted, label))

    # Incidence relations between r=1 and r=2 (if both levels present)
    levels_present = {r for r, _ in cell['levels']}
    if 1 in levels_present and 2 in levels_present:
        for rel, label in incidence_relations(n, 1):
            substituted = sp.expand(rel.subs(subs_map))
            results.append(('incidence', '1↔2', rel, substituted, label))

    return results


# ── Cell ideal via elimination ────────────────────────────────────────────────

def cell_ideal(cell, n, max_total_degree=4):
    """
    Compute the polynomial relations satisfied by the cell's Plücker coordinates.
    
    Strategy: introduce a symbolic p-variable for each (level, index-tuple) and
    set p - (polynomial in a's) = 0, then eliminate the a's via a Groebner basis.

    Returns a list of polynomials (relations) in the p-variables.
    """
    # Build the list of p-variables and their polynomial values in a's
    p_eqs = []
    p_vars = []
    a_vars = list(cell['tgt_free_vars'])

    for r, dets in cell['levels']:
        for idx_tuple, expr in dets:
            sym = plucker_symbol(r, tuple(i + 1 for i in idx_tuple))
            p_vars.append(sym)
            p_eqs.append(sym - expr)

    if not a_vars:
        # No free variables: each p_var is a constant, so substitution gives literal relations
        return [sp.simplify(eq) for eq in p_eqs if sp.simplify(eq) != 0]

    # Eliminate the a-variables to get relations purely among the p's
    try:
        gb = sp.groebner(p_eqs, *a_vars, *p_vars, order='lex')
        # Keep only polynomials in p_vars (with no a's)
        a_var_set = set(a_vars)
        relations = []
        for poly in gb:
            poly_expr = poly.as_expr()
            if not (poly_expr.free_symbols & a_var_set):
                if poly_expr != 0:
                    relations.append(poly_expr)
        return relations
    except Exception as e:
        return [f"Groebner failed: {e}"]


# ── Main ───────────────────────────────────────────────────────────────────────

num_str = input("Enter number string: ")

digits = [int(d) for d in num_str]
start, parts = 1, []
for d in digits:
    parts.append('(' + ''.join(str(i) for i in range(start, start + d)) + ')')
    start += d
print("Partition: " + ''.join(parts))

n_total  = sum(digits)
tableaux = young_tableaux_from_digits(num_str)
J        = jordan_from_digits(num_str)
valid    = generate_valid_fillings(tableaux, num_str)
cells    = fillings_to_schubert_cells(valid)

# Build cell data
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

# ── Section 1: universal Plücker relations ────────────────────────────────────

print("\n" + "=" * 80)
print("UNIVERSAL PLÜCKER RELATIONS for G(k, " + str(n_total) + ")")
print("=" * 80)

for k_level in range(2, n_total):
    rels = plucker_relations_in_level(n_total, k_level)
    if rels:
        print(f"\nLevel {k_level}  ({len(rels)} relations):")
        for rel, label in rels:
            print(f"  {rel} = 0   [{label}]")

inc_rels = incidence_relations(n_total, 1)
if inc_rels:
    print(f"\nIncidence relations (1 ↔ 2)  ({len(inc_rels)} relations):")
    for rel, label in inc_rels:
        print(f"  {rel} = 0   [{label}]")

# ── Section 2: verification per cell ──────────────────────────────────────────

print("\n" + "=" * 80)
print("VERIFICATION  —  do our computed Plücker coordinates satisfy the relations?")
print("=" * 80)

for cell in cell_data:
    print(f"\nC_{cell['perm']}")
    results = verify_plucker_relations(cell, n_total)

    all_pass     = all(r[3] == 0 for r in results)
    non_trivials = [r for r in results if r[3] == 0]   # all should be 0 if cell is on Grassmannian
    failures     = [r for r in results if r[3] != 0]

    if all_pass:
        print(f"  All {len(results)} relations satisfied (evaluate to 0).")
    if failures:
        print(f"  ✗ {len(failures)} relations FAILED:")
        for kind, lvl, rel, subbed, label in failures:
            print(f"    [{kind} {lvl}]  {rel} = 0   →  evaluates to {subbed}")

# ── Section 3: cell ideal computation ─────────────────────────────────────────

print("\n" + "=" * 80)
print("CELL IDEALS  —  polynomial relations among the Plücker coordinates of each cell")
print("=" * 80)

for cell in cell_data:
    print(f"\nC_{cell['perm']}  (parameter space A^{len(cell['tgt_free_vars'])})")
    ideal = cell_ideal(cell, n_total)
    if not ideal:
        print("  (no extra relations beyond universal Plücker)")
        continue
    for poly in ideal:
        print(f"  {poly} = 0")