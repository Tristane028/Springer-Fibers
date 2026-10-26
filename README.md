<h1>Springer Fibers: Algebraic Geometry &amp; Combinatorics</h1>

<p><strong>Authors:</strong> Tristan Endo, Alexander Ryan<br>
<strong>Advisor:</strong> Raymond O. Chou<br>
<strong>Institution:</strong> University of California, San Diego<br>
<strong>Year:</strong> 2026 to present</p>

<blockquote>
  ⚠️ Heads up: this repository is actively being updated as the research moves forward, so results, implementations, and documentation are all still subject to change.
</blockquote>

**[→ Explore the Springer fiber cell closures interactively](https://tristane028.github.io/Springer-Fibers/Visuals.html)**

---

## Overview

This project studies the geometry and combinatorics of Springer fibers in type A. For a nilpotent operator in Jordan form we enumerate the cells of its Springer fiber (the fiber intersected with each Schubert cell), compute their Plücker coordinates, and search for the **closure relations** between cells: which cells the closure of a given cell meets. The goal is to understand the irreducible components and how their closures are put together.

The repository holds four things:

- **Closure data** for 99 Jordan types (two, three and four rows, up to n = 10; a 100th, (2,3,5), is still running), with an explicit witness curve for every relation.
- **An interactive viewer** for exploring that data: the closure diagram in 3D, nested intersections, a searchable relation matrix, and "halving routes" that walk a closure down to a single cell.
- **The computation script**, rewritten to be fast enough for (3,3,3) and the larger three- and four-row cases.
- **Working notes** developing the ideal-theoretic method for computing an intersection exactly.

---

## Conventions

These changed in October 2026. Every data file now records `"labeling": "column-first"`, and data made under the old conventions is not mixed in.

**Partition key.** The digits are the row lengths of the Young diagram from top to bottom: `25` has a row of 2 boxes over a row of 5 (dashes separate parts once n ≥ 10, as in `1-9`).

**Labeled tableau.** Boxes are numbered down each column, leftmost column first:

```
33  →  [1][3][5]        321  →  [1][4][6]        2111  →  [1][5]
       [2][4][6]               [2][5]                      [2]
                               [3]                         [3]
                                                           [4]
```

**Nilpotent operator.** Each box sends its basis vector to the box on its left; the first column goes to 0. For 22 this is e₃ ↦ e₁, e₄ ↦ e₂.

**Blocks and cells.** The blocks are the rows of the labeled tableau, read left to right: (135)(246) for 33, (146)(25)(3) for 321, (15)(2)(3)(4) for 2111. A cell is indexed by a permutation in which every block appears in increasing order. Its dimension is the number of free parameters left after the span conditions X·Vⱼ ⊆ Vⱼ are imposed.

---

## Computed data

Files are grouped by number of rows and by n:

```
Two-Row-Springer-Fibers/N = 7/data_25.json
Three-Row-Springer-Fibers/N = 9/data_333.json
Four-Row-Springer-Fibers/N = 7/data_1222.json
```

Each `data_<partition>.json` contains:

| Field | Meaning |
|---|---|
| `partition`, `labeling` | the Jordan type and the labeling convention (`column-first`) |
| `cells[]` | `perm` (the permutation), `dim`, `free` (free parameters), `plucker` (Plücker coordinates at every level of the flag) |
| `edges[]` | one closure relation each: `src`, `tgt`, the limit variable `limit`, and the substitutions `srcSubs`, `tgtSubs` that witness it |
| `_complete` | present once every pair has been checked; an unfinished run instead keeps a `_checked` list so it can resume |

**99 partitions · 25,588 cells · 709,116 closure relations.**

### What a relation means

An edge `src → tgt` means the closure of the cell C_src **meets** the cell C_tgt. The witness is a one-parameter family of points in C_src, given by the substitutions with the limit variable sent to infinity, whose limit flag lies in C_tgt. Each source parameter is set to 0, ±1, ±T, ±T² or ±T³, and the target point uses 0, ±1 or ±T. For example, in (4,4):

```
C_12468357 → C_12354768
  source:  a₅ = a₃²,  a₆ = a₃
  target:  a₁ = −1,  a₂ = −1
  limit:   a₃ → ∞    (the Plücker vectors converge level by level, up to scale)
```

Two things to keep in mind:

- **Relations do not compose.** If the closure of A meets B and the closure of B meets C, the closure of A need not meet C, because "meets" includes partial intersections. The viewer keeps relations reached only through chains of curves separate and marks them as unproven.
- **The data is complete for what it found, not for everything.** Every recorded relation has a witness, but the search only tries the values above, so a relation that needs, say, T⁴ or a coefficient of 2 is missed. See the checks below.

### Two-row partitions

Components are iterated ℙ¹-bundles (Fung), all smooth.

| Partition | n | Blocks (rows of the labeled tableau) | Cells | Relations | Top dim |
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

*44 partitions · 2,024 cells · 13,339 relations.*


### Three-row partitions

| Partition | n | Blocks (rows of the labeled tableau) | Cells | Relations | Top dim |
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
| 2,2,5 | 9 | (1,4) (2,5) (3,6,7,8,9) | 756 | 17,487 | 6 |
| 2,3,4 | 9 | (1,4) (2,5,7) (3,6,8,9) | 1,260 | 51,323 | 7 |
| 3,3,3 | 9 | (1,4,7) (2,5,8) (3,6,9) | 1,680 | 103,378 | 9 |
| 4,3,2 | 9 | (1,4,7,9) (2,5,8) (3,6) | 1,260 | 49,039 | 7 |
| 5,2,2 | 9 | (1,4,7,8,9) (2,5) (3,6) | 756 | 17,401 | 6 |
| 5,3,1 | 9 | (1,4,6,8,9) (2,5,7) (3) | 504 | 8,212 | 5 |
| 6,2,1 | 9 | (1,4,6,7,8,9) (2,5) (3) | 252 | 2,263 | 4 |
| 7,1,1 | 9 | (1,4,5,6,7,8,9) (2) (3) | 72 | 290 | 3 |
| 1,1,8 | 10 | (1) (2) (3,4,5,6,7,8,9,10) | 90 | 370 | 3 |
| 1,2,7 | 10 | (1) (2,4) (3,5,6,7,8,9,10) | 360 | 3,401 | 4 |
| 1,3,6 | 10 | (1) (2,4,6) (3,5,7,8,9,10) | 840 | 15,173 | 5 |
| 1,4,5 | 10 | (1) (2,4,6,8) (3,5,7,9,10) | 1,260 | 35,772 | 6 |
| 2,2,6 | 10 | (1,4) (2,5) (3,6,7,8,9,10) | 1,260 | 31,578 | 6 |
| 5,3,2 | 10 | (1,4,7,9,10) (2,5,8) (3,6) | 2,520 | 117,494 | 7 |
| 5,4,1 | 10 | (1,4,6,8,10) (2,5,7,9) (3) | 1,260 | 34,734 | 6 |
| 6,2,2 | 10 | (1,4,7,8,9,10) (2,5) (3,6) | 1,260 | 31,468 | 6 |
| 6,3,1 | 10 | (1,4,6,8,9,10) (2,5,7) (3) | 840 | 14,967 | 5 |
| 7,2,1 | 10 | (1,4,6,7,8,9,10) (2,5) (3) | 360 | 3,370 | 4 |
| 8,1,1 | 10 | (1,4,5,6,7,8,9,10) (2) (3) | 90 | 371 | 3 |

*49 partitions · 21,824 cells · 625,565 relations.*

(2,3,5) is still being searched; its file will be added when the run finishes.


### Four-row partitions

| Partition | n | Blocks (rows of the labeled tableau) | Cells | Relations | Top dim |
|---|---:|---|---:|---:|---:|
| 1,1,1,2 | 5 | (1) (2) (3) (4,5) | 60 | 653 | 6 |
| 2,1,1,1 | 5 | (1,5) (2) (3) (4) | 60 | 671 | 6 |
| 1,1,2,2 | 6 | (1) (2) (3,5) (4,6) | 180 | 3,822 | 7 |
| 2,2,1,1 | 6 | (1,5) (2,6) (3) (4) | 180 | 3,675 | 7 |
| 1,2,2,2 | 7 | (1) (2,5) (3,6) (4,7) | 630 | 31,574 | 9 |
| 2,2,2,1 | 7 | (1,5) (2,6) (3,7) (4) | 630 | 29,817 | 9 |

*6 partitions · 1,740 cells · 70,212 relations.*

### Checks on the data

- **Bruhat order.** Each cell lies in a Schubert cell, so the closure of C_w can only meet cells C_v with v ≤ w in Bruhat order. The search skips every other pair, and no recorded relation violates this.
- **Counting cells in a closure.** Each cell contains exactly one torus-fixed point, so the number of cells a closure meets (counting itself) equals its Poincaré polynomial at q = 1.
  - Components of two-row Springer fibers are iterated ℙ¹-bundles (Fung), so a d-dimensional one meets exactly 2ᵈ cells.
  - In the two-row data **no closure meets more than 2^dim cells**. 1,884 of 2,024 cell closures, and 871 of 965 components, meet exactly 2^dim.
  - The rest come up short, which locates the relations the search missed: for example 26 of the 42 components of (5,5).
- **A certified comparison.** An exact intersection report for (3,2,2), computed from full and partial pieces with their ideals, agrees with every relation the search found. It adds further intersections in 68 of the 210 closures.
- **Which shapes have only smooth components.** By Fresse–Melnikov, every component is smooth exactly when the Jordan type is a hook, two rows, two rows plus one box, or (2,2,2). Every other shape has a singular component and, by Fresse, a component whose Poincaré polynomial is not palindromic. The viewer states which case each partition is in.
- **Earlier verification.** The old version of this README reported an exact match with the noncrossing-matching theorem of Goldwasser–Nadeem–Sun–Tymoczko for two-row partitions. That check was run under the previous labeling and operator and has not yet been redone for the column-first data.

---

## Halving routes (power-of-two degeneration)

Start from a cell's closure set S: the cell together with every cell its closure meets. One step picks a cell B in S and keeps S ∩ ({B} ∪ cells the closure of B meets). The aim is to reach a single cell while the population shrinks by a power of 2 at each step.

**The rule.** The divisors used are 128, 64, 32, 16, 8, 4 and 2:

- They are taken in decreasing order: a divisor may repeat or be skipped, but never goes back up (no 8 → 4 → 16).
- Each step lands on ⌊n/d⌋ or ⌈n/d⌉ cells.
- Dividing by d needs at least d cells.
- The route ends at one cell, which is always the last cell chosen.

**Result.** Every one of the 25,489 cells in the repository (99 complete Jordan types) with more than one cell in its closure has such a route, with rounding either way and even with rounding down only.

**Why powers of 2.** For two rows the reason is geometric. A component is a tower of ℙ¹-bundles, each layer doubling the number of fixed points. Fixing the top k layers at a fixed point gives a sub-tower with 2^(d−k) of them, so dividing by 2ᵏ is peeling off k layers. The data agrees: inside complete two-row closures, every one-step intersection has a power-of-two size.

**Where it gets harder.** With three or more rows the closure counts carry odd factors, coming from the column factors [c]! (for example [3]! = (1+q)(1+q+q²), giving 6 and 12). A closure of 6 cells shaped like the flag variety of ℂ³ has sub-intersections of sizes 1, 2, 4 and 6 but never 3.

- **Pure halving** (÷2 at every step, rounding either way) fails for many cells with three or more rows. For two rows it is nearly exact: 2 failing cells out of 1,980.
- **The 16 → 2 ladder** failed for closures of 80, 112 and 144 cells (16 × 5, 16 × 7, 16 × 9) in (4,3,2), (1,2,2,2) and (2,2,2,1).
- **Extending to 128** absorbs those odd factors through rounding. As n grows the top of the ladder will need to grow with the closures.

**The cube picture.** When a cell of dimension d has a closure meeting exactly 2^d cells, its half-size sub-closures S ∩ clo(B) give the closure coordinates in {0,1}^d: bit i of a cell is 1 when it lies outside the i-th half-size sub-closure, so the cell itself is 11…1. A halving step is then fixing one coordinate.

- In all 1,840 such two-row cells there are exactly d half-size sub-closures, and they give a bijection onto {0,1}^d.
- In 1,547 of them every sub-closure S ∩ clo(B) is a face of the cube ({c : code(c) ≤ code(B)}), so the closure lattice is the Boolean lattice.
- Cell dimension does not always equal cube level (the number of 1s): in many closures some cells sit at the "wrong" level, so the cube is not graded by dimension.
- A **tower route** halves at every step while the chosen cell drops exactly one dimension. Their number is at most d!, the number of chains of faces. Conjecture, not yet proved: it equals the number of linear extensions of a nesting forest, d!/∏(subtree sizes).

---

## Interactive viewer

[`Visuals.html`](https://tristane028.github.io/Springer-Fibers/Visuals.html) is a single page that reads the data files straight from this repository. It needs no build step and no manifest: a new `data_*.json` anywhere in the repo appears on the next page load.

- **Catalog.** Partitions grouped by number of rows and n, with cell counts and top dimension. Each shape sits directly above its reverse (rows read bottom to top): the lexicographically smaller one on top, e.g. for three rows at n = 9 the top line is (1,1,7) (1,2,6) (1,3,5) (2,2,5) (2,3,4) (3,3,3) and the bottom line (7,1,1) (6,2,1) (5,3,1) (5,2,2) (4,3,2). A shape equal to its reverse stands alone. When only one of a pair is in the repository, the other shows as a dashed “not computed yet” placeholder that cannot be opened (a link straight to it returns to the catalog); it becomes a normal entry as soon as its `data_*.json` is uploaded.
- **Jordan block card.** The labeled tableau, the arrows of the nilpotent operator, the blocks, and the matrix of X.
- **Diagram.** The closure diagram as a rotating 3D stack of discs, one per dimension, with a flat layered view as an alternative. Each dimension has its own color (up to dimension 11). "Related cells only" trims large partitions to the selected cell's neighbourhood.
- **Nested intersections.** Click a cell to see what its closure meets. Click one of those to look from inside both closures, and keep nesting. The chain lives in the URL, so a view can be shared.
- **Witnesses.** Each relation shows its curve, re-checked exactly in the browser with rational arithmetic.
- **Matrix.** Every ordered pair as a grid square: yes (curve found), probably (only through a chain), unresolved (below in Bruhat order with no curve), or no. Axis labels, a crosshair and a "does the closure of _ meet _?" box make a single relation quick to find.
- **Halving route animation.** From any cell, an animation of its 128 → … → 2 route: the chosen cell is ringed, its closure reaches out, and the cells that drop away vanish one by one. A counter and a log-scale "cells left" curve track the population. Routes come in a gradual style (halve unless forced) or with the fewest steps.
- **Coordinates tab.** For the selected cell: the nilpotent operator *X* (each box sent to the box on its left, boxes labeled column-first), the Schubert cell matrix *C* (pivot 1 in row w_j of column j, free entries above it, the relations from X·V_j ⊆ V_j already substituted), the product X·C, and every Plücker coordinate of every V_r = span of the first r columns, zeros included, written as one projective vector per V_r, [p_I : p_J : …] with the index sets in lexicographic order (zeros greyed, nonzeros highlighted, hover for the index, a Copy button for plain text); a labelled-table view with a switch to hide zeros is also available. *C* is read back from the stored Plücker vectors: with P the pivot rows of the first r−1 columns, entry C[s][r] = ± p_{P∪{s}}; this was checked entry by entry against `Computational Geometry.py`. Each matrix has a *Copy LaTeX* button; the Plücker table downloads as CSV for one cell or for every cell of the partition.
- **Routes tab.** Every cell's route as a table, with step counts, divisors and the chosen cells, plus a CSV download. Cells whose closure meets exactly 2^dim cells get a **⧉ cube** button (also reachable from the cell panel: *View closure as a cube*) that draws the closure as a d-cube, colour = cell dimension, dashed ring = dimension differs from its cube level. Two layouts: *by level* (rows = number of 1s in the code) and *box* (a square, a cube in perspective, a cube inside a cube for d = 4, two of those side by side for d = 5). A step control fades out each cut of the halving route. The table also counts tower routes and halving chains (a starred ⧉ marks a closure that is not a perfect cube).
- **Smoothness note.** For each partition, what Fresse–Melnikov and Fresse prove about its components.
- **Themes.** Light and dark, following the system setting.

---

## Repository contents

| Path | What it is |
|---|---|
| [`Visuals.html`](https://tristane028.github.io/Springer-Fibers/Visuals.html) | The interactive viewer |
| `Two-Row-Springer-Fibers/`, `Three-Row-Springer-Fibers/`, `Four-Row-Springer-Fibers/` | `data_<partition>.json`, one per Jordan type, in `N = n` subfolders |
| `Springer_Fiber_Closure (2).pdf`, `(3).pdf` | Working notes, described below |
| `Cell_Closures_for_Two-Row_Springer_Fibers.pdf` | Reference paper, see [References](#references) |

---

## Reproducing the computation

```
pip install numpy sympy python-flint
python3 "Computational Geometry.py" 333
```

The argument is the partition key (`333`, `2211`, `1-9`, …). The script:

1. builds every cell, imposing the span conditions and computing the Plücker coordinates at each level;
2. rules out every pair that Bruhat order forbids;
3. searches each remaining pair for a witness curve;
4. writes `data_<partition>.json` to your home folder.

Progress is saved every 30 seconds, so an interrupted run resumes where it stopped when you run the same command again. The finished file drops the resume bookkeeping and is saved as compact JSON.

**Settings** (environment variables):

| Variable | Default | Effect |
|---|---|---|
| `WORKERS` | all cores | number of worker processes; lower it on machines with little memory |
| `PRUNE_BRUHAT` | `1` | set to `0` to search pairs that Bruhat order rules out |
| `CHECKPOINT_SECS` | `30` | how often progress is saved |
| `EXACT_STRINGS` | `0` | `1` prints Plücker minors exactly as older versions did (slow for n ≥ 8) |
| `SIGN_LAZY`, `TERM_BUDGET`, `PROGRESS` | on | sign shortcut, memory cap per batch, per-pass progress lines |

**How it is made fast** (every change was checked to give exactly the same relations as before):

- **Cells.** Plücker minors are built level by level and simplified from their expanded form, with repeated minors cached: (3,3,2) builds in about 12 s instead of over 20 minutes, and (3,3,3) in about 85 s.
- **Level-by-level search.** Level r of the Plücker vector depends only on the first r columns. Parameters are therefore assigned in the order they first appear, and a partial substitution is dropped as soon as its levels so far match no remaining target.
- **Signs only when they matter.** A parameter's minus sign is tried only once it can change a leading coefficient, which saves about 3× on the largest cells.
- **Shared target index.** Every possible target point is computed once, hashed level by level, and shared between workers through memory-mapped files.
- **Earliest witness only.** For each target only the earliest witness is kept and then confirmed exactly, so memory stays flat. Batches are sized to cap memory per worker.

(3,3,3), with 1,680 cells and 294,838 pairs to search after Bruhat pruning, runs overnight on a laptop. The 42 top-dimensional cells take most of the time.

**Dependencies:** Python 3.9+, `numpy`, `sympy`, `python-flint`.

---

## Working notes

Two notes in this repository develop the ideal-theoretic method for computing a cell-closure intersection **exactly**, rather than certifying it pointwise. The computed relations say *whether* the closure of one cell meets another, but not *which subvariety* of the target the intersection is.

### Note I: saturation and the Rabinowitsch trick

📄 [View Note](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%282%29.pdf)

The note works through a small example: free variables *a*, *b*, *c* and the matching {(1,6),(2,3),(4,5)}. It sets up the core mechanism, saturation of an ideal *I* with respect to *f*, computed in practice via

$$(I : f^{\infty}) = (I + (1 - tf)) \cap R$$

so that the intersection is an elimination ideal of the augmented system. It ends with two questions: *which off-cell coordinates E<sub>I</sub> are needed in general, and what should one saturate by?*

### Note II: the N = 8 calculation

📄 [View Note](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Springer_Fiber_Closure%20%283%29.pdf)

The note is a full calculation for the two-row (4,4) fiber, with source matching M = {(1,8),(2,7),(3,4),(5,6)} and target F = {(1,2),(3,4),(5,6),(7,8)}. The result is

$$\overline{C_M} \cap C_F = \{s = p\} \subset C_F$$

which is a *proper subvariety* of the target cell, not the whole cell. The recipe:

1. Write the source cell in coordinates
2. Choose the target flag chart
3. Express the target-chart Plücker ratios on the source cell
4. Clear denominators and form the graph ideal
5. Saturate by the product of the denominators
6. Eliminate the source variables
7. **Only then** impose the equations cutting out the target cell

Order matters: imposing the target-cell equations too early is the main pitfall. Generically no finite point of the source cell lies in the target cell, so the intersection only appears once the closure is taken.

Saturation matters too. Eliminating from the unsaturated graph ideal gives only ⟨S − P⟩, which misses Y(P − Q) − X(P − R) = 0. Saturating by the pivot-minor product recovers ⟨S − P, Y(P − Q) − X(P − R)⟩. The note also shows that a few normal coordinates suffice: X = b/c and Y = b/d recover the full target-cell ideal. Macaulay2 code is included.

### Where this leads

These notes describe the exact computation behind the certified (3,2,2) report mentioned above. Making it routine would turn every relation in the viewer from "meets" into an explicit defining equation, and would fill in the relations the curve search misses.

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
