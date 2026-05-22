from sympy import symbols, Matrix, simplify
from sympy.parsing.sympy_parser import parse_expr
from itertools import combinations

#/usr/local/bin/python3.9 "/Users/christan065/springer-fibers-ray-chou/Plucker Coordinates.py"

def plucker_coordinates(matrix_entries, n, k, symbol_names=None):
    if symbol_names:
        syms = symbols(' '.join(symbol_names))
        if not isinstance(syms, tuple):
            syms = (syms,)
        local_dict = dict(zip(symbol_names, syms))
    else:
        local_dict = {}

    entries = [parse_expr(str(e), local_dict=local_dict) for e in matrix_entries]
    M = Matrix(n, k, entries)

    rows = list(range(n))
    all_levels = []

    for r in range(1, k):
        row_combos = list(combinations(rows, r))
        dets = []
        for rs in row_combos:
            # rows = subset, cols = first r columns
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
            idx = ''.join(str(r + 1) for r in rs)
            print(f"  p_{idx} = {d}")
            coords.append(d)
        print(f"  [{' : '.join(str(c) for c in coords)}]")

    print("\nPlücker coordinates:")
    parts = ['[' + ' : '.join(str(d) for _, d in dets) + ']'
             for _, dets in all_levels]
    print('(' + ', '.join(parts) + ')')


# --- Example: the 4x4 matrix from the problem ---
matrix_entries = [
    'x', '1', '0', '0', '0', '0',
    '0', '0', 'y', 'z', '1', '0',
    '0', '0', '0', 'y', '0', '1',
    '1', '0', '0', '0', '0', '0',
    '0', '0', '1', '0', '0', '0',
    '0', '0', '0', '1', '0', '0'
]

levels = plucker_coordinates(matrix_entries, n=6, k=6, symbol_names=['a', 'b'])
print_plucker(levels)
