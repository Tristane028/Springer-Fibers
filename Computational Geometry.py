import sys, json, time
import sympy as sp
from itertools import permutations, combinations, product as iproduct
import numpy as np
import os
from flint import fmpq, fmpq_poly

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Computational Geometry.py"

NUM_STR = sys.argv[1] if len(sys.argv) > 1 else '44'

SAVE_EVERY = 1   # checkpoint after this many computed pairs (1 = every pair)

# ── Build cells (UNCHANGED, SymPy) ────────────────────────────────────────────
def young_from(num_str): return [list(range(1, int(d)+1)) for d in num_str]
def jordan_from(num_str):
    digits = [int(d) for d in num_str]; n = sum(digits)
    J = np.zeros((n,n), dtype=int); idx = 0
    for k in digits:
        for i in range(k-1): J[idx+i, idx+i+1] = 1
        idx += k
    return sp.Matrix(J)
def blocks_of(num_str):
    digits = [int(d) for d in num_str]; bs, st = [], 1
    for d in digits: bs.append(list(range(st, st+d))); st += d
    return bs
def order_pres(perm, blocks):
    for b in blocks:
        if [perm.index(v) for v in b] != sorted([perm.index(v) for v in b]): return False
    return True
def valid_fillings(tab, num_str):
    n = sum(len(r) for r in tab); bs = blocks_of(num_str); v = []
    for p in permutations(range(1, n+1)):
        if not order_pres(p, bs): continue
        m, i = [], 0
        for row in tab: m.append(list(p[i:i+len(row)])); i += len(row)
        v.append(m)
    return v
def fillings_to_cells(valid):
    cs = []
    for mat in valid:
        flat = [v for row in mat for v in row]; n = len(flat)
        S = [[0]*n for _ in range(n)]
        for col, val in enumerate(flat): S[val-1][col] = 1
        cs.append(S)
    return cs
def fill_symbolic(cell):
    syms, idx = {}, 1; M = [r[:] for r in cell]
    for r in range(len(M)):
        for c in range(len(M)):
            if M[r][c] != 0: continue
            if not any(M[r][cc]==1 for cc in range(c)) and not any(M[rr][c]==1 for rr in range(r)):
                name = f"a{idx}"; M[r][c] = name; syms[name] = sp.symbols(name); idx += 1
    return M, syms
def span_check(C, XC):
    avars = sorted({s for s in C.free_symbols if s.name.startswith('a')},
                   key=lambda x: int(x.name[1:]))
    rel = {}
    for i in range(C.cols):
        Cs = C.subs(rel); XCs = XC.subs(rel)
        cols = Cs[:, :i+1]; tgt = XCs[:, i]
        alphas = sp.symbols(f"alpha0:{i+1}")
        sol = sp.solve(list(cols * sp.Matrix(alphas) - tgt), list(alphas) + avars, dict=True)
        if not sol: continue
        s = sol[0]
        for v in avars:
            if v in s and not any(a in s[v].free_symbols for a in alphas):
                rel[v] = s[v]
    return [v for v in avars if v not in rel], rel
def plucker(M, n, k):
    out = []
    for r in range(1, k):
        dets = []
        for rs in combinations(range(n), r):
            sub = M.extract(list(rs), list(range(r)))
            d = sub[0,0] if r==1 else sp.simplify(sub.det())
            dets.append((rs, d))
        out.append((r, dets))
    return out

# ── flint compilation: minors -> term lists for fast substitution ─────────────
def _to_fmpq(c):
    """sympy Integer/Rational -> flint fmpq."""
    r = sp.Rational(c)
    return fmpq(int(r.p), int(r.q))

def compile_levels(levels):
    """
    Pre-parse each Plücker minor into a flat term list so the inner loop never
    touches SymPy. Structure mirrors `levels`:
        compiled[level_index] = [ entry_terms, ... ]   (one per minor in the level)
        entry_terms          = [ (coeff: fmpq, {var_name: exponent}), ... ]
    """
    out = []
    for _, dets in levels:
        entries = []
        for _, d in dets:
            expr = sp.expand(sp.sympify(d))
            syms = sorted(expr.free_symbols, key=lambda s: int(s.name[1:]))
            terms = []
            if not syms:
                terms.append((_to_fmpq(expr), {}))
            else:
                poly = sp.Poly(expr, *syms)
                for monom, coeff in poly.terms():
                    exps = {syms[k].name: int(e) for k, e in enumerate(monom) if e}
                    terms.append((_to_fmpq(coeff), exps))
            entries.append(terms)
        out.append(entries)
    return out

# ── flint closure kernel ──────────────────────────────────────────────────────
_ONE = fmpq_poly([1])

def eval_entry(terms, subst):
    """Evaluate a compiled minor under subst (var_name -> fmpq_poly). Returns
    a single fmpq_poly in the limit variable."""
    acc = fmpq_poly([])
    for coeff, exps in terms:
        term = fmpq_poly([coeff])
        for var, e in exps.items():
            term = term * (subst[var] ** e)
        acc = acc + term
    return acc

def poly_lim(num, den):
    """lim_{x->oo} num/den as an fmpq, or None if it diverges (deg num > deg den).
    Canceling common factors never changes this, so no GCD is needed."""
    if num == 0:
        return fmpq(0)
    dn, dd = num.degree(), den.degree()
    if dn > dd: return None
    if dn < dd: return fmpq(0)
    return num[dn] / den[dd]

def try_step(src_c, tgt_c, src_subst, tgt_cache_entry):
    """
    Same logic as the SymPy try_step: per level, pick lambda from the first entry
    where source and target are both nonzero, require the scaled source AND the
    target to converge to the SAME finite value at every entry.

    tgt_cache_entry is the precomputed target side for this ts: a list (per level)
    of (tv_list, tl_list), or the sentinel None if the target diverges somewhere.
    """
    if tgt_cache_entry is None:
        return False
    for lvl_idx, sd_terms in enumerate(src_c):
        tv, tl = tgt_cache_entry[lvl_idx]
        sv = [eval_entry(t, src_subst) for t in sd_terms]
        # first index where both source and target are nonzero (lambda source)
        j = None
        for kk in range(len(sv)):
            if sv[kk] != 0 and tv[kk] != 0:
                j = kk; break
        for i in range(len(sv)):
            if j is None:
                sl = poly_lim(sv[i], _ONE)
            else:
                sl = poly_lim(tv[j] * sv[i], sv[j])   # (t_j / s_j) * s_i, as one ratio
            if sl is None or sl != tl[i]:
                return False
    return True

def clean(subs, levels):
    """Drop substitutions with no effect (key absent from the minors, or value ==
    key). UNCHANGED from the SymPy version; used only to format the recorded edge,
    so the JSON output stays byte-identical."""
    if not subs: return subs
    syms = set()
    for _, dets in levels:
        for _, d in dets: syms |= sp.sympify(d).free_symbols
    return {k:v for k,v in subs.items() if k in syms and sp.simplify(v-k)!=0}

def is_in_closure(src, tgt):
    """
    src, tgt are cell dicts carrying both 'compiled' (flint) and 'levels' (sympy,
    for clean/output). Returns (matched, limit_var, src_subs, tgt_subs) with the
    same semantics and output format as the SymPy version.
    """
    sf = src['free_vars']
    tf = tgt['tgt_free']
    if len(tf) > len(sf):
        return False, None, None, None

    src_c, tgt_c = src['compiled'], tgt['compiled']

    for lv in [None] + list(sf):
        other = [v for v in sf if v != lv]

        if lv is not None:
            T = fmpq_poly([0, 1])
            sym_sv = [sp.Integer(0), sp.Integer(1), sp.Integer(-1),
                      lv, lv**2, lv**3, -lv, -lv**2, -lv**3]
            fl_sv  = [fmpq_poly([0]), fmpq_poly([1]), fmpq_poly([-1]),
                      T, T**2, T**3, -T, -(T**2), -(T**3)]
            sym_tv = [sp.Integer(0), sp.Integer(1), sp.Integer(-1), lv, -lv]
            fl_tv  = [fmpq_poly([0]), fmpq_poly([1]), fmpq_poly([-1]), T, -T]
            base_subst = {lv.name: T}     # limit var stays live on the source side
        else:
            sym_sv = sym_tv = [sp.Integer(0), sp.Integer(1), sp.Integer(-1)]
            fl_sv  = fl_tv  = [fmpq_poly([0]), fmpq_poly([1]), fmpq_poly([-1])]
            base_subst = {}

        # Precompute the target side once per lv, keyed by the tuple of tgt indices.
        # (Target eval is independent of the source combo, so this is the analogue
        #  of the SymPy version's per-ts cache.)
        tgt_cache = {}
        for tc in iproduct(range(len(fl_tv)), repeat=len(tf)):
            ts_sym = {tf[m]: sym_tv[tc[m]] for m in range(len(tf))}
            if any(sp.simplify(v - k) == 0 for k, v in ts_sym.items()):
                tgt_cache[tc] = None          # identity sub -> skip (matches `continue`)
                continue
            tgt_subst = {tf[m].name: fl_tv[tc[m]] for m in range(len(tf))}
            per_level = []
            diverges = False
            for td_terms in tgt_c:
                tv = [eval_entry(t, tgt_subst) for t in td_terms]
                tl = [poly_lim(x, _ONE) for x in tv]
                if any(x is None for x in tl):   # whole target level diverges
                    diverges = True; break
                per_level.append((tv, tl))
            tgt_cache[tc] = None if diverges else per_level

        # Source combos (outer) x target combos (inner), original iteration order,
        # returning the first match so the recorded edge is identical.
        for sc in iproduct(range(len(fl_sv)), repeat=len(other)):
            src_subst = dict(base_subst)
            for m in range(len(other)):
                src_subst[other[m].name] = fl_sv[sc[m]]
            for tc in iproduct(range(len(fl_tv)), repeat=len(tf)):
                entry = tgt_cache[tc]
                if entry is None:
                    continue
                if try_step(src_c, tgt_c, src_subst, entry):
                    ss_sym = {other[m]: sym_sv[sc[m]] for m in range(len(other))}
                    ts_sym = {tf[m]: sym_tv[tc[m]] for m in range(len(tf))}
                    return (True, lv,
                            clean(ss_sym, src['levels']),
                            clean(ts_sym, tgt['levels']))
    return False, None, None, None

# ── Checkpoint I/O (UNCHANGED) ────────────────────────────────────────────────
def save_atomic(path, payload):
    """Write JSON to a temp file then rename over the target, so a crash mid-write
    cannot corrupt the existing checkpoint."""
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp, path)

# ── Run ───────────────────────────────────────────────────────────────────────
print(f"Computing for partition {NUM_STR}...", flush=True)
t0 = time.time()

J = jordan_from(NUM_STR)
fillings = valid_fillings(young_from(NUM_STR), NUM_STR)
raw_cells = fillings_to_cells(fillings)

cells_data = []
for i, cell in enumerate(raw_cells):
    M, syms = fill_symbolic(cell)
    C = sp.Matrix([[syms.get(v, v) for v in row] for row in M])
    XC = J * C
    for r in range(J.rows):
        if all(J[r,c]==0 for c in range(J.cols)):
            for c in range(XC.cols): XC[r,c] = 0
    free_vars, relations = span_check(C, XC)
    C_final = C.subs(relations)
    n, k = C_final.shape
    levels = plucker(C_final, n, k)
    perm = ''.join(str(v) for row in fillings[i] for v in row)
    psyms = set()
    for _, dets in levels:
        for _, d in dets: psyms |= sp.sympify(d).free_symbols
    active = [v for v in free_vars if v in psyms]
    cells_data.append({
        'perm':       perm,
        'levels':     levels,
        'compiled':   compile_levels(levels),     # flint term lists
        'free_vars':  free_vars,
        'tgt_free':   active,
        'plucker_strs': ['[' + ' : '.join(str(d) for _,d in dets) + ']' for _, dets in levels]
    })
    print(f"  cell {i+1}/{len(raw_cells)}: C_{perm}", flush=True)

# Build the JSON cell list once (cheap); edges accumulate below.
json_cells = [
    {
        'perm':    c['perm'],
        'dim':     len(c['tgt_free']),
        'free':    [str(v) for v in c['tgt_free']],
        'plucker': c['plucker_strs']
    }
    for c in cells_data
]

output_path = os.path.join(os.path.expanduser('~'), f'data_{NUM_STR}.json')

# ── Resume from a prior checkpoint, if present ─────────────────────────────────
edges   = []        # list of edge dicts
checked = set()     # "src>tgt" keys already evaluated

if os.path.exists(output_path):
    try:
        prior = json.load(open(output_path))
    except Exception:
        prior = None
    if prior and prior.get('partition') == NUM_STR:
        edges   = prior.get('edges', [])
        checked = set(prior.get('_checked', []))
        print(f"\nResuming from data_{NUM_STR}.json: "
              f"{len(checked)} pairs already evaluated, {len(edges)} edges found.",
              flush=True)
    else:
        print(f"\ndata_{NUM_STR}.json exists but is for a different partition; "
              f"starting fresh.", flush=True)

def build_out():
    return {
        'partition': NUM_STR,
        'cells':     json_cells,
        'edges':     edges,
        '_checked':  sorted(checked),
    }

# Write an initial checkpoint so the cells are on disk even before any edge.
save_atomic(output_path, build_out())

# ── Search closures, checkpointing as we go ────────────────────────────────────
print(f"\nSearching closures...", flush=True)
ordered_pairs = [
    (s, t)
    for s in cells_data
    for t in cells_data
    if s['perm'] != t['perm']
]
total_pairs = len(ordered_pairs)
since_save = 0

try:
    for done, (s, t) in enumerate(ordered_pairs, start=1):
        key = f"{s['perm']}>{t['perm']}"
        if key in checked:
            continue

        ok, lv, ss, ts = is_in_closure(s, t)
        checked.add(key)
        if ok:
            edges.append({
                'src': s['perm'],
                'tgt': t['perm'],
                'limit': str(lv) if lv else None,
                'srcSubs': {str(k): str(v) for k,v in (ss or {}).items()},
                'tgtSubs': {str(k): str(v) for k,v in (ts or {}).items()}
            })

        since_save += 1
        if since_save >= SAVE_EVERY:
            save_atomic(output_path, build_out())
            since_save = 0

        if done % 20 == 0 or ok:
            print(f"  {done}/{total_pairs} pairs checked, "
                  f"{len(edges)} edges ({time.time()-t0:.1f}s)", flush=True)

    save_atomic(output_path, build_out())
    print(f"\nWrote data_{NUM_STR}.json: {len(cells_data)} cells, "
          f"{len(edges)} edges in {time.time()-t0:.1f}s", flush=True)

except KeyboardInterrupt:
    save_atomic(output_path, build_out())
    print(f"\nInterrupted. Progress saved to data_{NUM_STR}.json "
          f"({len(checked)}/{total_pairs} pairs evaluated, {len(edges)} edges). "
          f"Re-run to resume.", flush=True)
    sys.exit(0)