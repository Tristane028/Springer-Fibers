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

## 1. Introduction

Springer fibers were introduced by Tonny Springer in 1976 as geometric objects whose cohomology carries representations of the symmetric group. Given a nilpotent matrix X, the Springer fiber is the set of complete flags V₁ ⊂ ⋯ ⊂ Vₙ with X·Vⱼ ⊆ Vⱼ for every j. Intersecting the fiber with Schubert cells breaks it into affine pieces, or cells, a result due to Spaltenstein and Shimomura. However, describing how these cells fit together is much harder. In particular, little is known about which cells the closure of a given cell meets.

For two Jordan blocks, Fung (2003) showed that every component is an iterated ℙ¹-bundle, and Goldwasser, Nadeem, Sun and Tymoczko (2025) described the cell closures completely using noncrossing matchings. Beyond two rows, no such description is known.

In this project, we compute the closure relations for 116 configurations with two, three and four rows, up to n = 11, which gives 35,614 cells and 1,144,312 relations. Furthermore, we verify the data against known results, find a rank condition that predicts the closures, and prove a cascade lemma that explains the structure of each cell. As a result, every cell can be written down directly from its permutation. We also build an interactive viewer for exploring the results.

---

## 2. Definitions

Consider a configuration, written as its row lengths from top to bottom, so `25` is a row of 2 boxes over a row of 5 (with dashes once a row has 10 or more boxes, as in `1-10`). We number the boxes down each column, starting with the leftmost column:

```
33  →  [1][3][5]        321  →  [1][4][6]
       [2][4][6]               [2][5]
                               [3]
```

Here, the nilpotent operator X sends each box to the box on its left, and the first column goes to 0. The rows of the tableau are the **groupings**. Each cell corresponds to a permutation in which every grouping appears in increasing order, and its dimension is the number of free variables left after we impose the flag.

---

## 3. Computed data

The data files live in `Two-`, `Three-` and `Four-Row-Springer-Fibers/N = n/`, and files above GitHub's upload limit are gzipped. Each file lists every cell (its permutation, dimension, free variables and Plücker coordinates) and every relation.

A relation `src → tgt` means the closure of C_src **meets** C_tgt. That is, there is a curve in C_src whose limit lands in C_tgt. For example, in (4,4):

```
C_12468357 → C_12354768
  source:  a₅ = a₃²,  a₆ = a₃
  target:  a₁ = −1,  a₂ = −1
  limit:   a₃ → ∞
```

Two remarks are in order. First, closure relations do not compose, so a relation reached only through a chain of curves is not proven. Second, the search only tries simple substitutions (0, ±1, ±T, ±T², ±T³), so some relations are missed. The rank test in Section 5 predicts most of these.

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

*54 configurations, 4,070 cells and 36,358 relations.*

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

*52 configurations, 25,604 cells and 774,938 relations.*

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

*10 configurations, 5,940 cells and 333,016 relations.*

---

## 4. Verification

- **Spanning check.** For every cell, each column of X·C lies in the span of the columns of C up to it, and every dimension matches the number of free variables.
- **Betti numbers.** The number of cells of each dimension matches the Betti numbers of the Springer fiber, computed independently from Kostka–Foulkes polynomials.
- **Bruhat order.** No relation goes against Bruhat order.
- **Explicit cells.** For all 35,614 cells, the explicit formula of Section 5.4 gives the same free variables and the same matrix as the spanning check.
- **Two rows.** No closure meets more than 2^dim cells, as Fung's theorem requires.
- **Exact comparison.** For (3,2,2), an exact intersection report computed with Gröbner bases agrees with every relation the search found.
- **Smoothness.** By Fresse and Melnikov, every component is smooth exactly when the configuration is a hook, two rows, two rows plus a box, or (2,2,2).

---

## 5. Results

### 5.1 The rank test

Define r_k(i, j) = dim(Xᵏ·V_j + V_i) − i. Since ranks can only drop in a limit,

> the closure of C_w can meet C_v only if v ≤ w in Bruhat order and every r_k(i, j) at v is at most its value on C_w.

None of the 1,144,312 relations breaks this condition. For two rows, it gives exactly 2^dim cells in all 4,070 closures, so it appears to be the full answer there. For (3,2,2), it matches the exact report in 209 of 210 closures, and the one exception is precisely the singular component. The predicted pairs are stored in `rank_<configuration>.json` and can be remade with `python3 rank_test.py`.

### 5.2 Cascading free variables

After we impose the flag, entries of C copy themselves from one column to the next column of the same grouping. This is because of the following lemma. If box q sits directly left of box p, in columns k′ and k, then

> X·C_k = C_k′ + (a combination of earlier columns whose 1 sits above row q).

In other words, the stack above the 1 in column k′ carries over into column k, unless an earlier arc boxes it in. Counting the variables that remain free gives the dimension of every cell,

> dim = Σ over boxes q of the smaller boxes after q and before the right neighbor of q,

which matches all 35,614 cells. A write-up with full proofs is under review.

### 5.3 The halving algorithm

Starting from a closure, we intersect with the closure of a chosen cell and divide the count by a power of 2, from 2¹⁰ = 1024 down to 2, until a single cell remains. Every one of the 35,498 cells with more than one cell in its closure has such a route. For two rows, each step peels off one ℙ¹ layer, so the closure behaves like a cube {0,1}^dim. With three or more rows, odd factors such as 3! = 6 appear, and rounding absorbs them.

### 5.4 Reading each cell off its permutation

Previously, the matrix of a cell came from solving the spanning check. Now it can be read directly from the permutation w, in the spirit of Section 4 of Goldwasser, Nadeem, Sun and Tymoczko. Index the entries by labels, so that E(r, q) is the entry in row r of the column whose 1 is in row q. Then each entry is one of the following.

- **Zero** if r comes before q in w.
- **Free** if q comes before r, and either q ends its row or r comes before R(q), the box to the right of q.
- **Copied** otherwise, from the entry one box to the right in both labels, E(R(r), R(q)), plus a correction.

The correction comes from the **blockers**, the labels b between r and q that come before R(q) in w. Each blocker adds y_b · D(b), where D(b) = E(R(b), R(q)) − E(b, q) is what slid in minus what is already there, and y_b is minus the entry above the 1 of b, after clearing it against the other blockers.

> E(r, q) = E(R(r), R(q)) + Σ over blockers b of y_b · D(b).

Nothing is solved, so each entry takes one pass. In other words, each free variable is copied diagonally along a **streak**, and a 1 in the way stops it and drags the entry above it along. Alexander Ryan found the same rule independently, as a rook that slides up and then left.

Furthermore, only blockers in a third grouping, different from those of r and q, can contribute. As such, **with two rows every entry is 0, 1 or a single variable**, which is the analogue of Theorem 4.10 of Goldwasser, Nadeem, Sun and Tymoczko. Polynomials first appear with three rows. In (3,3,3), 460 of the 2,520 dependent entries carry a correction.

To tell instantly whether an entry is a plain copy, list its blockers in a third grouping.

- If there are none, at this step and at every later slide, the entry is a copy. This is proven.
- If one has a nonzero entry above its 1 and a nonzero D(b), the entry is corrected. This never misses a correction, and its only false alarms are rare cancellations.

The formula is checked against every dependent entry in all 116 configurations. `flag_formula.py` prints the matrix of any cell, including configurations without data:

```
python3 flag_formula.py 333 362951487
```

### 5.5 Cutting a cell into smaller ones

Suppose the first m entries of w are 1, …, m in some order. Then the matrix is block-diagonal. The top block is a cell of the configuration formed by boxes 1, …, m, and the bottom block is a cell of the remaining boxes. Dimensions add, and the closure is the product of the two smaller closures. For example, (3,3,3) splits as 1 + 332, 11 + 322, 111 + 222, 211 + 221, and in reverse. In particular, the cell 326154879 is the (2,2,2) cell 326154 next to the (1,1,1) cell 213.

All 809 possible cuts in the data give the correct block form. Comparing closures across these cuts also turned up 64 relations that the search missed in eight smaller configurations, namely 321, 222, 322, 331, 133, 122, 114 and 112. Each one is forced by a curve found in a larger configuration, and the rank test had flagged all of them as allowed. They are listed in `missing_relations.json`.

---

## 6. Interactive viewer

[`Visuals.html`](https://tristane028.github.io/Springer-Fibers/Visuals.html) reads the data directly from this repository, so new files appear on the next page load.

- **Catalog.** Each configuration sits above its reverse, and missing ones appear as placeholders.
- **Diagram.** The closures in 3D, with nesting and an animation of the halving algorithm.
- **Matrix.** Every pair of cells, marked as found, reached through a chain, predicted by the rank test, or ruled out.
- **Routes.** The halving route of every cell, with a cube view for two-row closures.
- **Coordinates.** The operator X, the matrix before and after imposing the flag, the relations, the cascades and every Plücker coordinate, with variables named by position. The section *Reading C off the permutation* recomputes the matrix from w alone, as in Section 5.4. Each entry is marked free, copied along a streak, or corrected by blockers, and clicking an entry explains where it comes from.

---

## 7. Running the computation

```
pip install numpy sympy python-flint
python3 "Computational Geometry.py" 333
```

The script builds the cells, skips pairs that Bruhat order rules out, searches the rest, and writes `data_<configuration>.json` to the home folder. Progress is saved every 30 seconds, so rerunning the same command resumes an interrupted run. Versions of the script before November 15, 2026 could lose a relation when an interrupted run was resumed. This is now fixed, and computing a configuration again from scratch recovers anything that was lost. Two-row cases at n = 11 finish in under two minutes, while (3,3,3) runs overnight on a laptop.

| File | Purpose |
|---|---|
| `Computational Geometry.py` | the closure search |
| `rank_test.py` | generates the `rank_*.json` files |
| `flag_formula.py` | the matrix of any cell, read off its permutation |
| `missing_relations.json` | the 64 relations found by cutting cells |
| `build_manifest.py` | optional index of the data files |

---

## 8. Working notes

Two notes develop a method for computing an intersection **exactly** with ideals, rather than only deciding whether it is nonempty.

- **Note I: saturation and the Rabinowitsch trick.** 📄 [View](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%282%29.pdf) The note sets up saturation, (I : f^∞) = (I + (1 − tf)) ∩ R.
- **Note II: the N = 8 calculation.** 📄 [View](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%283%29.pdf) A full (4,4) example in which the intersection is a proper subvariety, {s = p}. In particular, the target-cell equations must be imposed only after elimination.

---

## Collaborators

- [Raymond Chou](https://github.com/raymondchou420)
- [Tristan Endo](https://github.com/Tristane028)
- [Alexander Ryan](https://github.com/Wave449)

## References

- Talia Goldwasser, Meera Nadeem, Garcia Sun, and Julianna Tymoczko. *Cell Closures for Two-Row Springer Fibers via Noncrossing Matchings.* [arXiv:2503.03941](https://arxiv.org/abs/2503.03941)
- Lucas Fresse and Anna Melnikov. *On the singularity of the irreducible components of a Springer fiber in sl(n).* [arXiv:0905.1617](https://arxiv.org/abs/0905.1617)
- Lucas Fresse. *Some characterizations of singular components of Springer fibers in the two-column case.* [arXiv:0909.4008](https://arxiv.org/abs/0909.4008)
- Francis Y. C. Fung. *On the topology of components of some Springer fibers and their relation to Kazhdan–Lusztig theory.* Advances in Mathematics 178 (2003). [arXiv:0020.4224](https://arxiv.org/pdf/math/0204224)
