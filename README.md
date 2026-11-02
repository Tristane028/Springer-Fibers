<h1>Springer Fibers: Algebraic Geometry &amp; Combinatorics</h1>

<p><strong>Authors:</strong> Tristan Endo, Alexander Ryan<br>
<strong>Advisor:</strong> Raymond O. Chou<br>
<strong>Institution:</strong> University of California, San Diego<br>
<strong>Year:</strong> 2026 to present</p>

<blockquote>
  ⚠️ This repository is actively being updated, so results and documentation may change.
</blockquote>

**[→ Explore the Springer fiber cell closures interactively](https://tristane028.github.io/Springer-Fibers/Visuals.html)**

---

## Overview

We study Springer fibers in type A. For each configuration we list the cells, compute their Plücker coordinates, and find which cells each closure meets. The goal is to understand the components and how their closures fit together.

**116 configurations · 35,614 cells · 1,144,312 closure relations**, up to n = 11.

---

## Conventions

- **Configuration key.** Row lengths from top to bottom: `25` is a row of 2 over a row of 5. Dashes once a row has 10+ boxes (`1-10`).
- **Labeling.** Boxes are numbered down each column, leftmost column first:

```
33  →  [1][3][5]        321  →  [1][4][6]
       [2][4][6]               [2][5]
                               [3]
```

- **Nilpotent operator.** X sends each box to the box on its left. The first column goes to 0.
- **Cells.** A cell is a permutation where each grouping (row of the tableau) stays in increasing order. Its dimension is the number of free variables left after imposing the flag (X·Vⱼ ⊆ Vⱼ).

---

## Data

Files live in `Two-`, `Three-` and `Four-Row-Springer-Fibers/N = n/`. Large files are gzipped (`.json.gz`).

| Field | Meaning |
|---|---|
| `cells[]` | permutation, dimension, free variables, Plücker coordinates |
| `edges[]` | one relation each, with the curve that proves it |
| `_complete` | set once every pair has been checked |

An edge `src → tgt` means the closure of C_src **meets** C_tgt. The proof is a curve in C_src whose limit lands in C_tgt. For example, in (4,4):

```
C_12468357 → C_12354768
  source:  a₅ = a₃²,  a₆ = a₃
  target:  a₁ = −1,  a₂ = −1
  limit:   a₃ → ∞
```

Relations don't compose, and the search only tries simple substitutions (0, ±1, ±T, ±T², ±T³), so a few relations get missed. The rank test below catches most of them.

### Two rows

| Configuration | n | Groupings | Cells | Relations | Top dim |
|---|---:|---|---:|---:|---:|
| 1,2 | 3 | (1) (2,3) | 3 | 2 | 1 |
| 2,1 | 3 | (1,3) (2) | 3 | 2 | 1 |
| 1,3 | 4 | (1) (2,3,4) | 4 | 3 | 1 |
| 2,2 | 4 | (1,3) (2,4) | 6 | 9 | 2 |
| 3,1 | 4 | (1,3,4) (2) | 4 | 3 | 1 |
| 1,4 | 5 | (1) (2,3,4,5) | 5 | 4 | 1 |
| 2,3 | 5 | (1,3) (2,4,5) | 10 | 19 | 2 |
| 3,2 | 5 | (1,3,5) (2,4) | 10 | 19 | 2 |
| 4,1 | 5 | (1,3,4,5) (2) | 5 | 4 | 1 |
| 1,5 | 6 | (1) (2,3,4,5,6) | 6 | 5 | 1 |
| 2,4 | 6 | (1,3) (2,4,5,6) | 15 | 32 | 2 |
| 3,3 | 6 | (1,3,5) (2,4,6) | 20 | 66 | 3 |
| 4,2 | 6 | (1,3,5,6) (2,4) | 15 | 32 | 2 |
| 5,1 | 6 | (1,3,4,5,6) (2) | 6 | 5 | 1 |
| 1,6 | 7 | (1) (2,3,4,5,6,7) | 7 | 6 | 1 |
| 2,5 | 7 | (1,3) (2,4,5,6,7) | 21 | 48 | 2 |
| 3,4 | 7 | (1,3,5) (2,4,6,7) | 35 | 145 | 3 |
| 4,3 | 7 | (1,3,5,7) (2,4,6) | 35 | 144 | 3 |
| 5,2 | 7 | (1,3,5,6,7) (2,4) | 21 | 48 | 2 |
| 6,1 | 7 | (1,3,4,5,6,7) (2) | 7 | 6 | 1 |
| 1,7 | 8 | (1) (2,3,4,5,6,7,8) | 8 | 7 | 1 |
| 2,6 | 8 | (1,3) (2,4,5,6,7,8) | 28 | 67 | 2 |
| 3,5 | 8 | (1,3,5) (2,4,6,7,8) | 56 | 262 | 3 |
| 4,4 | 8 | (1,3,5,7) (2,4,6,8) | 70 | 456 | 4 |
| 5,3 | 8 | (1,3,5,7,8) (2,4,6) | 56 | 261 | 3 |
| 6,2 | 8 | (1,3,5,6,7,8) (2,4) | 28 | 67 | 2 |
| 7,1 | 8 | (1,3,4,5,6,7,8) (2) | 8 | 7 | 1 |
| 1,8 | 9 | (1) (2,3,4,5,6,7,8,9) | 9 | 8 | 1 |
| 2,7 | 9 | (1,3) (2,4,5,6,7,8,9) | 36 | 89 | 2 |
| 3,6 | 9 | (1,3,5) (2,4,6,7,8,9) | 84 | 424 | 3 |
| 4,5 | 9 | (1,3,5,7) (2,4,6,8,9) | 126 | 1,033 | 4 |
| 5,4 | 9 | (1,3,5,7,9) (2,4,6,8) | 126 | 1,019 | 4 |
| 6,3 | 9 | (1,3,5,7,8,9) (2,4,6) | 84 | 423 | 3 |
| 7,2 | 9 | (1,3,5,6,7,8,9) (2,4) | 36 | 89 | 2 |
| 8,1 | 9 | (1,3,4,5,6,7,8,9) (2) | 9 | 8 | 1 |
| 1,9 | 10 | (1) (2,3,4,5,6,7,8,9,10) | 10 | 9 | 1 |
| 2,8 | 10 | (1,3) (2,4,5,6,7,8,9,10) | 45 | 114 | 2 |
| 3,7 | 10 | (1,3,5) (2,4,6,7,8,9,10) | 120 | 638 | 3 |
| 4,6 | 10 | (1,3,5,7) (2,4,6,8,9,10) | 210 | 1,963 | 4 |
| 5,5 | 10 | (1,3,5,7,9) (2,4,6,8,10) | 252 | 3,087 | 5 |
| 6,4 | 10 | (1,3,5,7,9,10) (2,4,6,8) | 210 | 1,946 | 4 |
| 7,3 | 10 | (1,3,5,7,8,9,10) (2,4,6) | 120 | 637 | 3 |
| 8,2 | 10 | (1,3,5,6,7,8,9,10) (2,4) | 45 | 114 | 2 |
| 9,1 | 10 | (1,3,4,5,6,7,8,9,10) (2) | 10 | 9 | 1 |
| 1,10 | 11 | (1) (2,3,4,5,6,7,8,9,10,11) | 11 | 10 | 1 |
| 2,9 | 11 | (1,3) (2,4,5,6,7,8,9,10,11) | 55 | 142 | 2 |
| 3,8 | 11 | (1,3,5) (2,4,6,7,8,9,10,11) | 165 | 911 | 3 |
| 4,7 | 11 | (1,3,5,7) (2,4,6,8,9,10,11) | 330 | 3,357 | 4 |
| 5,6 | 11 | (1,3,5,7,9) (2,4,6,8,10,11) | 462 | 7,168 | 5 |
| 6,5 | 11 | (1,3,5,7,9,11) (2,4,6,8,10) | 462 | 7,031 | 5 |
| 7,4 | 11 | (1,3,5,7,9,10,11) (2,4,6,8) | 330 | 3,338 | 4 |
| 8,3 | 11 | (1,3,5,7,8,9,10,11) (2,4,6) | 165 | 910 | 3 |
| 9,2 | 11 | (1,3,5,6,7,8,9,10,11) (2,4) | 55 | 142 | 2 |
| 10,1 | 11 | (1,3,4,5,6,7,8,9,10,11) (2) | 11 | 10 | 1 |

*54 configurations · 4,070 cells · 36,358 relations*

### Three rows

| Configuration | n | Groupings | Cells | Relations | Top dim |
|---|---:|---|---:|---:|---:|
| 1,1,2 | 4 | (1) (2) (3,4) | 12 | 34 | 3 |
| 2,1,1 | 4 | (1,4) (2) (3) | 12 | 35 | 3 |
| 1,2,2 | 5 | (1) (2,4) (3,5) | 30 | 163 | 4 |
| 2,2,1 | 5 | (1,4) (2,5) (3) | 30 | 159 | 4 |
| 1,1,4 | 6 | (1) (2) (3,4,5,6) | 30 | 106 | 3 |
| 1,2,3 | 6 | (1) (2,4) (3,5,6) | 60 | 413 | 4 |
| 2,2,2 | 6 | (1,4) (2,5) (3,6) | 90 | 1,118 | 6 |
| 3,2,1 | 6 | (1,4,6) (2,5) (3) | 60 | 413 | 4 |
| 4,1,1 | 6 | (1,4,5,6) (2) (3) | 30 | 107 | 3 |
| 1,1,5 | 7 | (1) (2) (3,4,5,6,7) | 42 | 157 | 3 |
| 1,2,4 | 7 | (1) (2,4) (3,5,6,7) | 105 | 823 | 4 |
| 1,3,3 | 7 | (1) (2,4,6) (3,5,7) | 140 | 1,593 | 5 |
| 2,2,3 | 7 | (1,4) (2,5) (3,6,7) | 210 | 3,623 | 6 |
| 3,2,2 | 7 | (1,4,7) (2,5) (3,6) | 210 | 3,570 | 6 |
| 3,3,1 | 7 | (1,4,6) (2,5,7) (3) | 140 | 1,559 | 5 |
| 4,2,1 | 7 | (1,4,6,7) (2,5) (3) | 105 | 822 | 4 |
| 5,1,1 | 7 | (1,4,5,6,7) (2) (3) | 42 | 158 | 3 |
| 1,1,6 | 8 | (1) (2) (3,4,5,6,7,8) | 56 | 218 | 3 |
| 1,2,5 | 8 | (1) (2,4) (3,5,6,7,8) | 168 | 1,433 | 4 |
| 1,3,4 | 8 | (1) (2,4,6) (3,5,7,8) | 280 | 4,038 | 5 |
| 2,2,4 | 8 | (1,4) (2,5) (3,6,7,8) | 420 | 8,660 | 6 |
| 2,3,3 | 8 | (1,4) (2,5,7) (3,6,8) | 560 | 17,010 | 7 |
| 3,3,2 | 8 | (1,4,7) (2,5,8) (3,6) | 560 | 16,159 | 7 |
| 4,2,2 | 8 | (1,4,7,8) (2,5) (3,6) | 420 | 8,600 | 6 |
| 4,3,1 | 8 | (1,4,6,8) (2,5,7) (3) | 280 | 3,962 | 5 |
| 5,2,1 | 8 | (1,4,6,7,8) (2,5) (3) | 168 | 1,426 | 4 |
| 6,1,1 | 8 | (1,4,5,6,7,8) (2) (3) | 56 | 219 | 3 |
| 1,1,7 | 9 | (1) (2) (3,4,5,6,7,8,9) | 72 | 289 | 3 |
| 1,2,6 | 9 | (1) (2,4) (3,5,6,7,8,9) | 252 | 2,280 | 4 |
| 1,3,5 | 9 | (1) (2,4,6) (3,5,7,8,9) | 504 | 8,327 | 5 |
| 1,4,4 | 9 | (1) (2,4,6,8) (3,5,7,9) | 630 | 14,010 | 6 |
| 2,2,5 | 9 | (1,4) (2,5) (3,6,7,8,9) | 756 | 17,487 | 6 |
| 2,3,4 | 9 | (1,4) (2,5,7) (3,6,8,9) | 1,260 | 51,323 | 7 |
| 3,3,3 | 9 | (1,4,7) (2,5,8) (3,6,9) | 1,680 | 103,378 | 9 |
| 4,3,2 | 9 | (1,4,7,9) (2,5,8) (3,6) | 1,260 | 49,039 | 7 |
| 4,4,1 | 9 | (1,4,6,8) (2,5,7,9) (3) | 630 | 13,764 | 6 |
| 5,2,2 | 9 | (1,4,7,8,9) (2,5) (3,6) | 756 | 17,401 | 6 |
| 5,3,1 | 9 | (1,4,6,8,9) (2,5,7) (3) | 504 | 8,212 | 5 |
| 6,2,1 | 9 | (1,4,6,7,8,9) (2,5) (3) | 252 | 2,263 | 4 |
| 7,1,1 | 9 | (1,4,5,6,7,8,9) (2) (3) | 72 | 290 | 3 |
| 1,1,8 | 10 | (1) (2) (3,4,5,6,7,8,9,10) | 90 | 370 | 3 |
| 1,2,7 | 10 | (1) (2,4) (3,5,6,7,8,9,10) | 360 | 3,401 | 4 |
| 1,3,6 | 10 | (1) (2,4,6) (3,5,7,8,9,10) | 840 | 15,173 | 5 |
| 1,4,5 | 10 | (1) (2,4,6,8) (3,5,7,9,10) | 1,260 | 35,772 | 6 |
| 2,2,6 | 10 | (1,4) (2,5) (3,6,7,8,9,10) | 1,260 | 31,578 | 6 |
| 2,3,5 | 10 | (1,4) (2,5,7) (3,6,8,9,10) | 2,520 | 121,599 | 7 |
| 5,3,2 | 10 | (1,4,7,9,10) (2,5,8) (3,6) | 2,520 | 117,494 | 7 |
| 5,4,1 | 10 | (1,4,6,8,10) (2,5,7,9) (3) | 1,260 | 34,734 | 6 |
| 6,2,2 | 10 | (1,4,7,8,9,10) (2,5) (3,6) | 1,260 | 31,468 | 6 |
| 6,3,1 | 10 | (1,4,6,8,9,10) (2,5,7) (3) | 840 | 14,967 | 5 |
| 7,2,1 | 10 | (1,4,6,7,8,9,10) (2,5) (3) | 360 | 3,370 | 4 |
| 8,1,1 | 10 | (1,4,5,6,7,8,9,10) (2) (3) | 90 | 371 | 3 |

*52 configurations · 25,604 cells · 774,938 relations*

### Four rows

| Configuration | n | Groupings | Cells | Relations | Top dim |
|---|---:|---|---:|---:|---:|
| 1,1,1,2 | 5 | (1) (2) (3) (4,5) | 60 | 653 | 6 |
| 2,1,1,1 | 5 | (1,5) (2) (3) (4) | 60 | 671 | 6 |
| 1,1,2,2 | 6 | (1) (2) (3,5) (4,6) | 180 | 3,822 | 7 |
| 2,2,1,1 | 6 | (1,5) (2,6) (3) (4) | 180 | 3,675 | 7 |
| 1,1,2,3 | 7 | (1) (2) (3,5) (4,6,7) | 420 | 11,652 | 7 |
| 1,2,2,2 | 7 | (1) (2,5) (3,6) (4,7) | 630 | 31,574 | 9 |
| 2,2,2,1 | 7 | (1,5) (2,6) (3,7) (4) | 630 | 29,817 | 9 |
| 3,2,1,1 | 7 | (1,5,7) (2,6) (3) (4) | 420 | 11,592 | 7 |
| 1,2,2,3 | 8 | (1) (2,5) (3,6) (4,7,8) | 1,680 | 122,975 | 9 |
| 3,2,2,1 | 8 | (1,5,8) (2,6) (3,7) (4) | 1,680 | 116,585 | 9 |

*10 configurations · 5,940 cells · 333,016 relations*

---

## Checks

- **Spanning check.** Every cell satisfies X·Vⱼ ⊆ Vⱼ, and every dimension matches its free variables.
- **Betti numbers.** Cell counts by dimension match the Kostka–Foulkes formula for every configuration.
- **Bruhat order.** No relation goes against Bruhat order.
- **Two rows.** No closure meets more than 2^dim cells, as Fung's theorem requires.
- **(3,2,2).** An exact intersection report agrees with every relation found.
- **Smoothness.** By Fresse–Melnikov, all components are smooth exactly for hooks, two rows, two rows plus a box, and (2,2,2).

---

## Results

### Rank test

With r_k(i, j) = dim(Xᵏ·V_j + V_i) − i:

> the closure of C_w can meet C_v only if v ≤ w in Bruhat order and every r_k(i, j) at v is at most its value on C_w.

- None of the 1,144,312 relations breaks it.
- For two rows it gives exactly 2^dim cells in all 4,070 closures, so it looks like the full answer there.
- For (3,2,2) it matches the exact report in 209 of 210 closures. The one miss is the singular component.

`rank_<configuration>.json` lists the pairs it allows that the search hasn't proven. Remake them with `python3 rank_test.py`.

### Cascading free variables

If box q sits left of box p, in columns k′ and k:

> X·C_k = C_k′ + (earlier columns with their 1 above row q).

So the stack above the 1 in column k′ carries over into column k unless an earlier arc boxes it in. Counting what stays free gives the dimension of every cell:

> dim = Σ over boxes q of the smaller boxes after q and before q's right neighbor.

This matches all 35,614 cells. A write-up with proofs is in progress.

### Halving algorithm

Pick a cell in the closure, intersect with its closure, and divide the count by a power of 2 (2¹⁰ = 1024 down to 2, rounding either way) until one cell is left.

- Every one of the 35,498 cells with more than one cell in its closure has a route.
- For two rows each step peels off one ℙ¹ layer, so the closure looks like a cube {0,1}^dim.
- With three or more rows, factors like 3! = 6 appear, and rounding absorbs them.

---

## Viewer

[`Visuals.html`](https://tristane028.github.io/Springer-Fibers/Visuals.html) reads the data straight from this repository. New files show up on the next page load.

- **Catalog.** Each configuration sits above its reverse. Missing ones show as placeholders.
- **Diagram.** The closures in 3D, with nesting and the halving algorithm animation.
- **Matrix.** Every pair at a glance: found, through a chain, predicted by the rank test, or no.
- **Routes.** Every cell's halving route, plus a cube view for two-row closures.
- **Coordinates.** X, the matrix before and after imposing the flag, the relations, the cascades and every Plücker coordinate. Variables are named by position, a_{row,col}.

---

## Running the computation

```
pip install numpy sympy python-flint
python3 "Computational Geometry.py" 333
```

It builds the cells, skips pairs Bruhat order rules out, searches the rest, and writes `data_<configuration>.json` to your home folder. Progress saves every 30 seconds, so rerunning the same command resumes. Set `WORKERS` to limit cores.

Two-row cases at n = 11 take under two minutes. (3,3,3) runs overnight on a laptop.

| File | What it is |
|---|---|
| `Computational Geometry.py` | the closure search |
| `rank_test.py` | makes the `rank_*.json` files |
| `build_manifest.py` | optional index of the data files |

---

## Working notes

Two notes work out how to compute an intersection **exactly** with ideals, not just whether it's nonempty.

- **Note I: saturation and the Rabinowitsch trick.** 📄 [View](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%282%29.pdf) Sets up saturation, (I : f^∞) = (I + (1 − tf)) ∩ R.
- **Note II: the N = 8 calculation.** 📄 [View](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%283%29.pdf) A full (4,4) example where the intersection is a proper subvariety, {s = p}. Key lesson: impose the target-cell equations only after eliminating.

---

## Collaborators

- [Raymond Chou](https://github.com/raymondchou420)
- [Tristan Endo](https://github.com/Tristane028)
- [Alexander Ryan](https://github.com/Wave449)

## References

- Talia Goldwasser, Meera Nadeem, Garcia Sun, and Julianna Tymoczko. *Cell Closures for Two-Row Springer Fibers via Noncrossing Matchings.* [arXiv:2503.03941](https://arxiv.org/abs/2503.03941) · 📄 [View Paper](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Cell_Closures_for_Two-Row_Springer_Fibers.pdf)
- Lucas Fresse and Anna Melnikov. *On the singularity of the irreducible components of a Springer fiber in sl(n).* [arXiv:0905.1617](https://arxiv.org/abs/0905.1617)
- Lucas Fresse. *Some characterizations of singular components of Springer fibers in the two-column case.* [arXiv:0909.4008](https://arxiv.org/abs/0909.4008)
- Francis Y. C. Fung. *On the topology of components of some Springer fibers and their relation to Kazhdan–Lusztig theory.* Advances in Mathematics 178 (2003).
