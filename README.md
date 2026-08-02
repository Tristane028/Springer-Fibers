<h1>Springer Fibers — Algebraic Geometry &amp; Combinatorics</h1>

<p><strong>Authors:</strong> Tristan Endo, Alexander Ryan<br>
<strong>Advisor:</strong> Raymond O. Chou<br>
<strong>Institution:</strong> University of California, San Diego<br>
<strong>Year:</strong> 2026 — Present</p>

<blockquote>
  ⚠️ This repository is actively being updated as new research is conducted. Results, implementations, and documentation are subject to change.
</blockquote>

**[→ Explore the Springer fiber cell complex interactively](https://tristane028.github.io/Springer-Fibers/Visuals.html)**

---

## Overview

This project investigates the geometry and combinatorics of Springer fibers. We examine the structure of flag varieties under nilpotent operators defined by Jordan blocks, compute valid Young tableau fillings, and analyze the corresponding Schubert cells. The goal is to characterize the irreducible components of Springer fibers through Plücker coordinates, Springer span checks, closure relations, and explicit defining equations.

The repository currently holds three things: **computed closure data** for 14 partitions, an **interactive visualizer** for exploring it, and **working notes** developing the ideal-theoretic method for computing intersections exactly.

---

## Computed data

Each `data_<partition>.json` records every Schubert cell of that Springer fiber — its permutation label, dimension, free parameters, and Plücker coordinates at every level of the flag — together with the directed closure relations between cells and an explicit witness for each one.

The partition key gives the Jordan block sizes in order, and the blocks partition {1, …, N} into consecutive runs. So `25` means two Jordan blocks of sizes 2 and 5 acting on ℂ⁷, with blocks {(1,2), (3,4,5,6,7)} and Young diagram

```
▢ ▢ ▢ ▢ ▢
▢ ▢
```

**14 partitions · 352 cells · 1,662 closure relations.**

### Two-row partitions

These are the cases covered by the noncrossing-matching theory, where each cell corresponds to a standard noncrossing matching and the closure order is combinatorially predictable.

| Partition | Shape | Blocks | Cells | Relations | Top dim |
|---|---|---|---:|---:|---:|
| 1,2 | ▢▢ / ▢ | { (1), (2,3) } | 3 | 2 | 1 |
| 1,3 | ▢▢▢ / ▢ | { (1), (2,3,4) } | 4 | 3 | 1 |
| 1,5 | ▢▢▢▢▢ / ▢ | { (1), (2,3,4,5,6) } | 6 | 5 | 1 |
| 2,2 | ▢▢ / ▢▢ | { (1,2), (3,4) } | 6 | 9 | 2 |
| 2,4 | ▢▢▢▢ / ▢▢ | { (1,2), (3,4,5,6) } | 15 | 32 | 2 |
| 2,5 | ▢▢▢▢▢ / ▢▢ | { (1,2), (3,4,5,6,7) } | 21 | 48 | 2 |
| 3,2 | ▢▢▢ / ▢▢ | { (1,2,3), (4,5) } | 10 | 19 | 2 |
| 3,3 | ▢▢▢ / ▢▢▢ | { (1,2,3), (4,5,6) } | 20 | 67 | 3 |
| 3,4 | ▢▢▢▢ / ▢▢▢ | { (1,2,3), (4,5,6,7) } | 35 | 146 | 3 |
| 4,4 | ▢▢▢▢ / ▢▢▢▢ | { (1,2,3,4), (5,6,7,8) } | 70 | 473 | 4 |

*10 partitions · 190 cells · 804 relations.*

### Three-row partitions

Three Jordan blocks. These lie outside the noncrossing-matching theory, so there is no independent combinatorial check on their closure relations — but they reach higher dimension relative to their size, which is part of what makes them interesting.

| Partition | Shape | Blocks | Cells | Relations | Top dim |
|---|---|---|---:|---:|---:|
| 2,2,1 | ▢▢ / ▢▢ / ▢ | { (1,2), (3,4), (5) } | 30 | 172 | 5 |
| 3,2,1 | ▢▢▢ / ▢▢ / ▢ | { (1,2,3), (4,5), (6) } | 60 | 423 | 5 |
| 4,1,1 | ▢▢▢▢ / ▢ / ▢ | { (1,2,3,4), (5), (6) } | 30 | 106 | 3 |
| 5,1,1 | ▢▢▢▢▢ / ▢ / ▢ | { (1,2,3,4,5), (6), (7) } | 42 | 157 | 3 |

*4 partitions · 162 cells · 858 relations.*

Note that (2,2,1) and (3,2,1) reach dimension 5 on only 5 and 6 boxes, while the two-row (4,4) needs 8 boxes to reach dimension 4 — the three-row fibers are considerably more intricate per box.

### Witnesses

Every edge carries the data that realizes it: the substitution applied to the source parameters, the values fixed on the target, and which variable is sent to infinity. For example, in (4,4):

```
C_56781234 → C_15672348
  source:  a₄ = a₆⁴,  a₇ = a₆²,  a₈ = a₆³
  target:  a₃ = 0,  a₆ = 0,  a₉ = −1
  limit:   a₆ → ∞
```

### Verification against the noncrossing-matching theorem

For two-row partitions, Goldwasser–Nadeem–Sun–Tymoczko (Theorem 6.9) give the closure as a union of *cuts* of the source's noncrossing matching: cutting an arc swaps the letters of the associated {B,T}-word at the arc's endpoints. This gives an independent combinatorial prediction of the closure order, which we use as a regression test on the computed geometry.

**All ten two-row partitions agree exactly** — every predicted relation present, no spurious extras:

| Partition | Computed | Predicted | |
|---|---:|---:|:--|
| 1,2 · 1,3 · 1,5 | 2 · 3 · 5 | 2 · 3 · 5 | ✓ |
| 2,2 · 2,4 · 2,5 | 9 · 32 · 48 | 9 · 32 · 48 | ✓ |
| 3,2 · 3,3 · 3,4 | 19 · 67 · 146 | 19 · 67 · 146 | ✓ |
| 4,4 | 473 | 473 | ✓ |

The cell dimensions match too: dimension equals the number of arcs in the corresponding matching, confirmed for every cell of every two-row partition.

The three-row partitions lie outside the theorem's scope and have no comparable combinatorial check; their completeness rests on the search having run over all ordered cell pairs, which the recorded bookkeeping confirms it did.

#### A correction this check produced

The cross-check initially reported (4,4) at 471 relations against 473 predicted, with the shortfall isolated to a single dimension transition (4 → 3). Both missing relations turned out to originate from `C_56781234` — the fully nested matching `BBBBTTTT`, the top-dimensional cell — and both require a **degree-4** substitution (`a₆⁴`, `a₄⁴`). The search grid at the time capped polynomial substitutions at degree 3, so no witness could be found. Widening the grid recovered both immediately; each was then verified independently by symbolic limits across all seven Plücker levels. Both are now included.

The lesson generalizes: if degree 3 was insufficient at (4,4), it is likely insufficient for larger partitions, and the substitution degree should be treated as a parameter to check rather than a fixed constant.

---

## Interactive visualizer

`Visuals.html` is a self-contained page — all data embedded, no build step, no dependencies — for exploring the cell complex.

- **Partition picker.** Grouped into two-row and three-row families, so the two regimes stay visually distinct.
- **Block structure.** Each partition shows its Young diagram with the numbered boxes, alongside the block decomposition of {1, …, N} — for (2,5), the diagram above together with { (1,2), (3,4,5,6,7) }.
- **Overview.** All cells of a partition as a 3-D web, coloured by dimension, height encoding dimension, with the closed point at the bottom. Relations are drawn undirected here so the global structure reads clearly.
- **Cell view.** Open any cell to see its parameter space 𝔸<sup>d</sup>, its Plücker coordinates level by level, and its outgoing closure relations as arrows.
- **Path view.** Select a relation to see just that pair, framed with source left and target right, and watch the degeneration animate as a stream of particles drawn into the destination. The substitutions and the limit variable stay pinned on screen after it completes.

Notation renders in LaTeX style throughout — italic variables with true subscripts and superscripts, so `p^(2)_12 = a₆²` appears as it would in a paper rather than as ASCII.

---

## Working notes

Two notes in this repository develop the ideal-theoretic method for computing a cell-closure intersection **exactly**, rather than certifying it pointwise. This is the natural next step beyond the current data: the computed relations answer *whether* one cell meets another's closure, but not *which subvariety* of the target that intersection is.

### `Springer_Fiber_Closure__2_.pdf` — saturation and the Rabinowitsch trick

Works a small example (free variables *a*, *b*, *c*; matching {(1,6),(2,3),(4,5)}) and establishes the mechanism: saturation of an ideal *I* with respect to *f*, computed in practice via

$$(I : f^{\infty}) = (I + (1 - tf)) \cap R$$

so that the intersection is obtained as an elimination ideal from the augmented system. It closes with the two questions that drive the second note: *which off-cell coordinates E<sub>I</sub> are needed in general, and what should one saturate by in general?*

### `Springer_Fiber_Closure__3_.pdf` — the N = 8 calculation

A full worked calculation for the two-row (4,4) fiber, with source matching M = {(1,8),(2,7),(3,4),(5,6)} and target F = {(1,2),(3,4),(5,6),(7,8)}. The result is

$$\overline{C_M} \cap C_F = \{s = p\} \subset C_F$$

— a *proper subvariety* of the target cell, not the whole cell. The note gives a seven-step recipe:

1. Write the source cell in coordinates
2. Choose the target flag chart
3. Express the target-chart Plücker ratios on the source cell
4. Clear denominators and form the graph ideal
5. Saturate by the product of the denominators
6. Eliminate the source variables
7. **Only then** impose the equations cutting out the target cell

Order matters: imposing the target-cell equations too early is the central pitfall, since generically no finite point of the source cell lies in the target cell — the intersection appears only after taking the closure.

Saturation is likewise not cosmetic. Eliminating from the unsaturated graph ideal yields only ⟨S − P⟩, missing the relation Y(P − Q) − X(P − R) = 0; saturating by the pivot-minor product recovers the correct chart closure ⟨S − P, Y(P − Q) − X(P − R)⟩. The note also shows that a small set of normal coordinates suffices in practice — here X = b/c and Y = b/d recover the full target-cell ideal — and includes Macaulay2 code for the computation.

### Where this leads

The visualizer currently treats each closure relation as binary: reachable or not, with one witness point. The notes show the real content is a defining equation. Surfacing that equation per edge would turn the graph from a reachability diagram into a geometric one, and is the clearest direction for extending both the computation and the interface.

---

## Reproducing the computation

The closure search takes a number string (e.g. `"221"`) where each digit is a Jordan block size, generates all order-preserving tableau fillings, converts each to a Schubert cell matrix, runs Springer span checks to find free variables, computes Plücker coordinates, and then tests every ordered pair of cells for a closure relation by searching over polynomial substitutions with a limit variable sent to infinity, checking projective convergence of the Plücker vectors.

Output is written incrementally to `data_<partition>.json` and the run is resumable.

**Performance notes.** The pair search dominates the runtime and is worth attention on larger partitions:

- Target-side evaluation depends only on the target cell and the limit variable, never on the source — so it should be memoized across sources rather than rebuilt inside every pair.
- Within a source substitution, the source Plücker vector is constant across all target combinations; evaluating it once per level instead of once per target combination is the single largest saving.
- A zero-pattern prefilter rules out most pairs before any symbolic work: wherever the source's Plücker coordinate vanishes identically, the target's must be able to vanish too, and a nonzero *constant* never can. This is a necessary condition, so it never rejects a true relation — on (2,2,1) it eliminates about 79% of pairs.
- Checkpointing after every pair is quadratic I/O on large partitions, since the whole file is rewritten each time. Batch it.

**Dependencies:** Python 3.9+, `numpy`, `sympy`, and `python-flint` for the compiled polynomial arithmetic.

```
pip install numpy sympy python-flint
```

---

## Collaborators

- [Raymond Chou](https://github.com/raymondchou420)
- [Tristan Endo](https://github.com/Tristane028)
- [Alexander Ryan](https://github.com/Wave449)

## References

Talia Goldwasser, Meera Nadeem, Garcia Sun, and Julianna Tymoczko. *Cell Closures for Two-Row Springer Fibers via Noncrossing Matchings.* [arXiv:2503.03941](https://arxiv.org/abs/2503.03941) · [📄 View local copy](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Cell_Closures_for_Two-Row_Springer_Fibers.pdf)
