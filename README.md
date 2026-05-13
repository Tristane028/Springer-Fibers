# Springer Fibers — Algebraic Geometry & Combinatorics

**Authors:** Tristan Endo, Alexander Ryan  
**Advisor:** Raymond O. Chou  
**Institution:** University of California, San Diego  
**Year:** 2026 — Present

> ⚠️ This repository is actively being updated as new research is conducted. Results, implementations, and documentation are subject to change.

## Overview

This project investigates the geometry and combinatorics of Springer fibers. We examine the structure of flag varieties under nilpotent operators defined by Jordan blocks, compute valid Young tableau fillings, and analyze the corresponding Schubert cells. The goal of this research is to characterize the irreducible components of Springer fibers through Plücker coordinates and Springer span checks.

## What the Code Does

The script takes a number string as input (e.g. `"221"`) where each digit defines the size of a Jordan block. From there it does the following:

- Constructs the Young tableau and Jordan matrix from the input
- Generates all valid fillings of the tableau that respect the order-preserving condition within each Jordan block
- Converts each valid filling into a Schubert cell matrix with symbolic entries
- Runs Springer span checks to determine which entries are free variables and which are constrained
- Computes Plücker coordinates for each cell using minors of the matrix
- Outputs a full summary table of all cells with their tableau, matrix, free variables, relations, and Plücker coordinates

## Running the Script

```bash
python "Combo Research.py"
```

The script will prompt you to enter a number string:

```
Enter number string: 221
```

## Dependencies

- Python 3.9
- `numpy`
- `sympy`

Install with:

```bash
pip install numpy sympy
```

## Collaborators

- [Tristan Endo](https://github.com/Tristane028)
- [Alexander Ryan](https://github.com/Wave449)

## References

Talia Goldwasser, Meera Nadeem, Garcia Sun, and Julianna Tymoczko. *Cell Closures for Two-Row Springer Fibers via Noncrossing Matchings.* 📄 [View Paper](https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Cell_Closures_for_Two-Row_Springer_Fibers.pdf)
