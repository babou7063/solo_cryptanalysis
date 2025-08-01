# Cryptanalysis Project

## Sources and Inspirations

This project implements several variants of the Pollard's rho algorithm for solving the discrete logarithm problem on elliptic curves. It includes basic walks, additive walks, and an optimized version over the field \( \mathbb{Z}/(2^{128} - 3) \mathbb{Z} \).

### Theoretical References

- D. J. Bernstein, T. Lange, P. Schwabe. "On the correct use of the negation map in the Pollard rho method". 2010. https://dl.acm.org/doi/10.5555/1964658.1964669
- Christophe Petit, "Cryptanalysis", Cours INFO-F514, Université libre de Bruxelles, 2025.
- Algorithms from classical cryptanalysis: Pollard's rho, additive walk, arithmetic on elliptic curves.
- [Pollard's rho algorithm for logarithms – Wikipedia](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm_for_logarithms)
- [The Hare and the Hedgehog – Wikipedia](https://en.wikipedia.org/wiki/The_Hare_and_the_Hedgehog)
- [Handbook of Applied Cryptography – Chapter 3: Number-Theoretic Reference Problems](https://cacr.uwaterloo.ca/hac/about/chap3.pdf) (Menezes, van Oorschot, Vanstone)

### Public Repositories Consulted

- [github.com/StackeredSAS/Pollard_Rho](https://github.com/StackeredSAS/Pollard_Rho/tree/main) – Used as a point of reference for structural comparison and implementation on small groups.

### Tools Used

- [Windsurf.ai – Generate Docstring](https://windso.rs/) was used to generate initial docstring templates, which were then reviewed and adapted manually.
- [OpenAI Codex](https://openai.com/blog/openai-codex) was used to generate or assist in the drafting of basic function scaffolding, especially for helper functions and parsing logic.





############################################
# Cryptanalysis Project: Pollard's Rho Algorithms for Elliptic Curve Discrete Logarithm

This project implements several variants of Pollard's rho algorithm for solving the discrete logarithm problem (DLP) on elliptic curves. The implementation includes basic walks, additive walks, and an optimized version over the finite field ℤ/(2^128 - 3)ℤ with negation maps. Given points P and Q on an elliptic curve, where Q = kP for some unknown scalar k, the goal is to efficiently find k.


## Core Components

### Elliptic Curve Infrastructure
- **`elliptic_curve.py`**: Standard elliptic curve operations over finite fields ℤ/pℤ
- **`mod128_minus3.py`**: Specialized arithmetic for the field ℤ/(2^128 - 3)ℤ with 10-coefficient representation

### Algorithm Implementations
1. **Basic Rho** (`basic_rho.py`) : Classical Floyd's cycle detection
2. **Additive Walk** (`additive_walk_rho.py`) : Precomputed table-based walks
3. **Negation Map Optimization** (`rho_mod128.py`) : Advanced variant with fruitless cycle detection


### Prerequisites
- Python 3.8+
- matplotlib (for visualizations)
- Standard library modules: random, hashlib, time

### Basic Usage

```python
...
```

## Algorithm Comparison

| Algorithm | Time Complexity | Space Complexity | Special Features |
|-----------|----------------|------------------|------------------|
| Basic Rho | O(√n) | O(1) | Simple Floyd cycle detection |
| Additive Walk | O(√n) | O(r) | Precomputed table, better parallelization |
| Negation Map | O(√n/2) | O(r) | Fruitless cycle detection, canonical forms |

Where n is the order of the base point and r is the precomputed table size.

# TODO -> rajouter ma comparaison