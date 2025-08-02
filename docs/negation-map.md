# Negation Map Pollard's Rho Algorithm

## Overview

The Negation Map variant of Pollard's rho algorithm enhances the additive walk approach by using canonical point representations to handle the negation map efficiently. This method reduces fruitless cycles by mapping each point P to its canonical form |P| (choosing between P and -P), making the walk more robust and achieving better practical performance on elliptic curve discrete logarithm problems.

## Algorithm Description

### Core Innovation

The negation map approach combines the additive walk strategy with canonical point representations. Each point is mapped to its canonical form using y-coordinate parity:

```
|P| = P if y is even, else |P| = -P
Wi+1 = |Wi + Rh(Wi)|  where h(Wi) is a hash function
Rj = cj·P + dj·Q      (precomputed with random cj, dj)
```

The algorithm maintains coefficient tracking:
```
ai+1 = ai + ch(Wi) (mod order)
bi+1 = bi + dh(Wi) (mod order)
```

### Key Features

1. **Canonical Form Mapping**: Uses y-coordinate parity to choose between P and -P
2. **Fruitless Cycle Detection**: Identifies and escapes 2-cycles and 4-cycles
3. **Distinguished Points**: Uses trailing zero bits in x-coordinate for collision detection
4. **Additive Walk**: Precomputed table of random linear combinations for better randomization

### Fruitless Cycle Handling

The algorithm detects common problematic patterns:
- **2-cycles**: Wi-1 = Wi-3 (detected every ~2√r steps)
- **4-cycles**: Wi-1 = Wi-5 
- **Escape mechanism**: Doubles the minimum point by x-coordinate and resets coefficients


## Complexity Analysis

### Time Complexity
- **Expected**: O(√n) where n is the order of point P
- **Theoretical bound**: ~√(πn/4) group operations (improved by negation map)

### Space Complexity  
- **O(r)** where r is the precomputed table size (default: 32-64)
- **Distinguished points storage**: O(d) where d is the number of distinguished points

### Performance Characteristics
- **Table size r**: Typically 32-64 for small curves, larger for production use
- **Distinguished point density**: 1/8 (3 zero bits) provides good balance
- **Cycle detection frequency**: Every 2√r steps to minimize overhead

## Usage Examples

### Basic Example
```python
from elliptic_curve import *
from rho_negative_map import *

# Setup curve and points
curve = EllipticCurve(a=2, b=3, p=101)
P = Point(1, 2, curve)
order = curve.find_order(P)

k_secret = 7
Q = curve.scalar_mul(k_secret, P)

# Solve using negation map rho
solver = SmallECRho(curve, P, Q, order, r=32)
k_found = solver.solve_ecdlp(max_walks=500)
print(f"Found k = {k_found}")
```