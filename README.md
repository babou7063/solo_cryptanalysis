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

