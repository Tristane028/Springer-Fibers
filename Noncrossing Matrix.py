import re
import sys
import sympy as sp

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Noncrossing Matrix.py"

# ── labels ────────────────────────────────────────────────────────────────────
def arc_label(i):
    """1 -> 'a', 2 -> 'b', ... 26 -> 'z', 27 -> 'a1', 28 -> 'a2', ... (rare)."""
    if 1 <= i <= 26:
        return chr(ord('a') + i - 1)
    return f"a{i}"          # fallback for very large matchings


# ── validation ──────────────────────────────────────────────────────────────--
def validate_matching(pairs):
    """
    Normalise and validate. Returns arcs sorted by left endpoint: list of
    (l, r) with l < r. Raises ValueError describing the first problem found.
    """
    arcs = []
    for p in pairs:
        if len(p) != 2:
            raise ValueError(f"each pair needs exactly 2 endpoints, got {p!r}")
        a, b = int(p[0]), int(p[1])
        if a == b:
            raise ValueError(f"an arc cannot connect a point to itself: {p!r}")
        arcs.append((min(a, b), max(a, b)))

    endpoints = [x for (l, r) in arcs for x in (l, r)]
    n = len(endpoints)
    if sorted(endpoints) != list(range(1, n + 1)):
        raise ValueError(
            f"endpoints must be exactly 1..{n} each used once; got {sorted(endpoints)}")

    # crossing test: arcs (l1,r1),(l2,r2) cross iff l1<l2<r1<r2 (strictly interleaved)
    s = sorted(arcs)
    for i in range(len(s)):
        l1, r1 = s[i]
        for j in range(i + 1, len(s)):
            l2, r2 = s[j]
            if l1 < l2 < r1 < r2:
                raise ValueError(
                    f"arcs {(l1, r1)} and {(l2, r2)} cross "
                    f"(need nested or disjoint, not interleaved)")
    return s


# ── construction ──────────────────────────────────────────────────────────────
def build(pairs, as_sympy=True):
    """
    Build the matrix for a noncrossing matching.

    Returns a dict with:
      'arcs'   : arcs sorted by left endpoint, index i (0-based) is arc i+1
      'labels' : {arc_index_1based: label}
      'top'    : m x n list of lists (entries are 0, 1, or a label string)
      'bottom' : m x n list of lists
      'grid'   : full n x n list of lists (top stacked over bottom)
      'matrix' : sympy.Matrix of the full grid (if as_sympy) with labels as Symbols
    """
    arcs = validate_matching(pairs)
    m = len(arcs)
    n = 2 * m
    labels = {i + 1: arc_label(i + 1) for i in range(m)}

    # 1-based l, r and ancestor lookups
    L = {i + 1: arcs[i][0] for i in range(m)}
    R = {i + 1: arcs[i][1] for i in range(m)}

    def ancestors(i):
        """arcs strictly containing arc i, ordered innermost -> outermost."""
        anc = [j for j in range(1, m + 1)
               if j != i and L[j] < L[i] and R[j] > R[i]]
        anc.sort(key=lambda j: L[j], reverse=True)   # largest l = innermost first
        return anc

    top    = [[0] * n for _ in range(m)]
    bottom = [[0] * n for _ in range(m)]
    current_row = {}     # arc -> row (1-based) it has most recently descended to

    # process events left to right; at any column at most one start or one end
    start_at = {L[i]: i for i in range(1, m + 1)}
    end_at   = {R[i]: i for i in range(1, m + 1)}

    for col in range(1, n + 1):
        if col in start_at:
            i = start_at[col]
            chain = [i] + ancestors(i)          # top -> bottom
            length = len(chain)
            bottom_row = i                       # outermost lands in row i
            top_row = bottom_row - (length - 1)  # arc i's own row
            for k, arc in enumerate(chain):
                row = top_row + k
                top[row - 1][col - 1] = labels[arc]
                current_row[arc] = row
            # bottom block: arc i's pivot at its start column
            bottom[i - 1][col - 1] = 1
        elif col in end_at:
            i = end_at[col]
            top[current_row[i] - 1][col - 1] = 1

    grid = [row[:] for row in top] + [row[:] for row in bottom]

    out = {'arcs': arcs, 'labels': labels, 'top': top, 'bottom': bottom, 'grid': grid}

    if as_sympy:
        syms = {lab: sp.Symbol(lab) for lab in labels.values()}
        M = sp.zeros(n, n)
        for r in range(n):
            for c in range(n):
                v = grid[r][c]
                M[r, c] = syms[v] if isinstance(v, str) else sp.Integer(v)
        out['matrix'] = M
    return out


# ── pretty printing ─────────────────────────────────────────────────────────--
def render(pairs):
    res = build(pairs, as_sympy=False)
    top, bottom = res['top'], res['bottom']
    n = len(top[0])
    w = max(len(str(x)) for row in res['grid'] for x in row)
    def fmt_row(row):
        return "[" + " ".join(str(x).rjust(w) for x in row) + "]"
    lines = []
    for row in top:
        lines.append(fmt_row(row))
    lines.append("-" * len(lines[0]))
    for row in bottom:
        lines.append(fmt_row(row))
    return "\n".join(lines)


def render_arcs(pairs, ascii=False):
    """
    Draw the matching as nested labelled arcs above the number line, e.g.

        ┌────────(a)────────┐
        │   ┌(b)┐   ┌(c)┐   │
        │   │   │   │   │   │
        1   2   3   4   5   6

    Outermost arcs sit on top; each arc is drawn at a row equal to its nesting
    depth, with its label centred on the bar. Set ascii=True for .-\\,_,| glyphs.
    """
    res = build(pairs, as_sympy=False)
    arcs, labels = res['arcs'], res['labels']
    m = len(arcs); n = 2 * m
    Lp = {i + 1: arcs[i][0] for i in range(m)}
    Rp = {i + 1: arcs[i][1] for i in range(m)}

    TL, TR, HR, VR = ('.', '.', '_', '|') if ascii else ('┌', '┐', '─', '│')

    def depth(i):
        return sum(1 for j in range(1, m + 1)
                   if j != i and Lp[j] < Lp[i] and Rp[j] > Rp[i])
    dep = {i: depth(i) for i in range(1, m + 1)}
    D = max(dep.values(), default=0)

    labtxt = {i: f"({labels[i]})" for i in range(1, m + 1)}
    maxlab = max((len(t) for t in labtxt.values()), default=3)
    s = max(4, maxlab + 1)
    col = lambda i: (i - 1) * s
    width = (n - 1) * s + 1
    tick = D + 1
    g = [[' '] * width for _ in range(tick + 1)]

    for i in range(1, m + 1):
        d, lc, rc = dep[i], col(Lp[i]), col(Rp[i])
        g[d][lc], g[d][rc] = TL, TR
        for c in range(lc + 1, rc):
            g[d][c] = HR
        lab = labtxt[i]
        start = lc + 1 + ((rc - lc - 1) - len(lab)) // 2
        for k, ch in enumerate(lab):
            g[d][start + k] = ch
        for rr in range(d + 1, tick):
            if g[rr][lc] == ' ': g[rr][lc] = VR
            if g[rr][rc] == ' ': g[rr][rc] = VR
    for i in range(1, n + 1):
        g[tick][col(i)] = VR

    lines = [''.join(r).rstrip() for r in g]
    numline = [' '] * width
    for i in range(1, n + 1):
        for k, ch in enumerate(str(i)):
            if col(i) + k < width:
                numline[col(i) + k] = ch
    lines.append(''.join(numline).rstrip())
    return '\n'.join(lines)


def show(pairs, ascii=False):
    """
    Arc diagram and matrix on one shared grid, so point i, number i, and matrix
    column i all sit in the same character column. Returns the combined string.
    """
    res = build(pairs, as_sympy=False)
    arcs, labels, grid = res['arcs'], res['labels'], res['grid']
    m = len(arcs); n = 2 * m
    Lp = {i + 1: arcs[i][0] for i in range(m)}
    Rp = {i + 1: arcs[i][1] for i in range(m)}
    TL, TR, HR, VR = ('.', '.', '_', '|') if ascii else ('┌', '┐', '─', '│')

    def depth(i):
        return sum(1 for j in range(1, m + 1)
                   if j != i and Lp[j] < Lp[i] and Rp[j] > Rp[i])
    dep = {i: depth(i) for i in range(1, m + 1)}
    D = max(dep.values(), default=0)

    labtxt = {i: f"({labels[i]})" for i in range(1, m + 1)}
    maxlab = max((len(t) for t in labtxt.values()), default=3)
    maxw = max((len(str(x)) for row in grid for x in row), default=1)
    s = max(4, maxlab + 1, maxw + 2)
    LM = 1                       # left margin reserves a column for '['
    P = lambda i: LM + (i - 1) * s
    W = P(n) + maxw + 1          # room for the closing ']'

    # arc diagram
    tick = D + 1
    g = [[' '] * W for _ in range(tick + 1)]
    for i in range(1, m + 1):
        d, lc, rc = dep[i], P(Lp[i]), P(Rp[i])
        g[d][lc], g[d][rc] = TL, TR
        for c in range(lc + 1, rc):
            g[d][c] = HR
        lab = labtxt[i]
        start = lc + 1 + ((rc - lc - 1) - len(lab)) // 2
        for k, ch in enumerate(lab):
            g[d][start + k] = ch
        for rr in range(d + 1, tick):
            if g[rr][lc] == ' ': g[rr][lc] = VR
            if g[rr][rc] == ' ': g[rr][rc] = VR
    for i in range(1, n + 1):
        g[tick][P(i)] = VR
    lines = [''.join(r).rstrip() for r in g]

    num = [' '] * W
    for i in range(1, n + 1):
        for k, ch in enumerate(str(i)):
            num[P(i) + k] = ch
    lines.append(''.join(num).rstrip())
    lines.append('')

    # matrix, columns aligned to the same P(i)
    def matrow(row):
        cells = [' '] * W
        cells[0] = '['
        for i in range(1, n + 1):
            for k, ch in enumerate(str(row[i - 1])):
                cells[P(i) + k] = ch
        cells[P(n) + maxw] = ']'
        return ''.join(cells)
    for r in grid[:m]:
        lines.append(matrow(r))
    lines.append('-' * (P(n) + maxw + 1))
    for r in grid[m:]:
        lines.append(matrow(r))
    return '\n'.join(lines)


# ── BT-string <-> pairs ───────────────────────────────────────────────────────
def parse_pairs(line):
    """Pull every integer out of the line and pair them in order, so
    '(1,8),(2,7),(3,6),(4,5)' -> [(1,8),(2,7),(3,6),(4,5)]. Kept for library use."""
    nums = [int(x) for x in re.findall(r'\d+', line)]
    if not nums:
        raise ValueError("no positions found — type pairs like (1,8),(2,7),(3,6),(4,5).")
    if len(nums) % 2 != 0:
        raise ValueError("odd number of endpoints — every position needs a partner.")
    return [(nums[i], nums[i + 1]) for i in range(0, len(nums), 2)]


def parse_bt_string(s):
    """
    Turn a BT-string (B = begin arc, T = terminate arc) into matching pairs.
    Scans left to right keeping a stack of open B positions; each T closes the
    most recent open B. A balanced BT-string is always a noncrossing matching.

        'BBTBTT' -> [(2, 3), (4, 5), (1, 6)]   (positions are the string indices)
    """
    s = re.sub(r'\s+', '', s).upper()
    if not s:
        raise ValueError("empty string — type something like BBTBTT.")
    bad = sorted({ch for ch in s if ch not in 'BT'})
    if bad:
        raise ValueError(f"only B and T are allowed; found {bad}.")

    stack, pairs = [], []
    for idx, ch in enumerate(s, start=1):
        if ch == 'B':
            stack.append(idx)
        else:  # T
            if not stack:
                raise ValueError(f"the T at position {idx} has no open B before it.")
            pairs.append((stack.pop(), idx))
    if stack:
        raise ValueError(f"{len(stack)} B(s) never closed (open at position(s) {stack}).")
    return pairs


def pairs_to_bt_string(pairs):
    """Inverse: noncrossing matching -> BT-string (B at every left endpoint)."""
    arcs = validate_matching(pairs)
    n = 2 * len(arcs)
    lefts = {l for (l, r) in arcs}
    return ''.join('B' if p in lefts else 'T' for p in range(1, n + 1))


# ── interactive entry ─────────────────────────────────────────────────────────
def ask_total_positions():
    while True:
        try:
            raw = input("Total Positions: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            sys.exit(0)
        try:
            n = int(raw)
        except ValueError:
            print("  Enter a whole even number, e.g. 8.\n")
            continue
        if n < 2 or n % 2 != 0:
            print("  Total positions must be even and at least 2 (4, 6, 8, ...).\n")
            continue
        return n


def ask_bt_string(n):
    while True:
        try:
            line = input("Noncrossing BT-String: ")
        except (EOFError, KeyboardInterrupt):
            print()
            sys.exit(0)
        cleaned = re.sub(r'\s+', '', line).upper()
        try:
            if len(cleaned) != n:
                raise ValueError(
                    f"string must have exactly {n} characters (one per position); "
                    f"got {len(cleaned)}.")
            pairs = parse_bt_string(cleaned)
            build(pairs)            # full validation (belt and suspenders)
        except ValueError as e:
            print(f"  {e}\n  Try again.\n")
            continue
        return pairs


def run_interactive():
    n = ask_total_positions()
    pairs = ask_bt_string(n)
    print()
    print(show(pairs))
    res = build(pairs, as_sympy=False)
    legend = "   ".join(f"{res['labels'][i + 1]}=({a},{b})"
                         for i, (a, b) in enumerate(res['arcs']))
    print("\n" + legend)
    return pairs


if __name__ == "__main__":
    run_interactive()