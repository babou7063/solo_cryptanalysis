# Negation Map Pollard's Rho Algorithm

## Overview

The Negation Map variant of Pollard’s rho exploits the curve symmetry (x,y)↦(x,−y) to reduce the search space by a factor of 2. In theory, this yields a √2 speedup but in practice constant factors (branching, canonicalization, cycle handling) can cancel the gain.

Our implementation enforces a canonical form (selecting between \(P\) and \(-P\) depending on the parity of y) and adds detection of fruitless cycles. On small curves in Python, results confirm that the theoretical advantage is masked: the negation-map walk is not faster than the additive walk, although it remains robust.


## Algorithm Description

### Core idea

Each step maps the current point to its canonical representative before updating:

```
|P| = P if y is even, else -P  
Wi+1 = |Wi + Rh(Wi)|  
```

with:
- `Rj = cj·P + dj·Q` precomputed in a jump table,  
- coefficients `(ai, bi)` updated modulo the subgroup order.  

### Features

1. **Canonical form mapping** – ensures consistency across equivalent points  
2. **Cycle handling** – detects 2- and 4-cycles; escapes by perturbing the walk (e.g. doubling a recent point)  
3. **Distinguished points** – collision detection based on x-coordinate bits  
4. **Precomputed additive walk** – amortizes scalar multiplications


## Why not mod 2^{128}−3 ?

Early attempts targeted F_2^128−3, but this was abandoned because:

- The curve order ∣E(Fp)∣ and subgroup order n are unknown making modular reduction and DP rates infeasible.  
- Computing ∣E(Fp)∣ with SEA has complexity O~(p1/4), out of reach for p≈2^128.  
- Python arithmetic adds prohibitive overhead at this scale.  

Instead, experiments were restricted to small primes where subgroup orders are known, allowing reliable comparisons of variants.


## Usage Example

```
from elliptic_curve import *  
from rho_negative_map import *  
  
curve = EllipticCurve(a=2, b=3, p=101)  
P = Point(1, 2, curve)  
order = curve.find_order(P)  
  
k_secret = 7  
Q = curve.scalar_mul(k_secret, P)  
  
solver = SmallECRho(curve, P, Q, order, r=32)  
k_found = solver.solve_ecdlp(max_walks=500)  
print(f"Found k = {k_found}")  
```