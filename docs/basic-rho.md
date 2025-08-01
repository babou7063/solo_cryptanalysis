# Basic Pollard's Rho Algorithm

## Overview

The basic Pollard's rho algorithm is the foundational method for solving the discrete logarithm problem on elliptic curves. This implementation uses Floyd's cycle detection algorithm (tortoise and hare) to find collisions in a pseudorandom walk on the elliptic curve group.

## Algorithm Description

### Core Principle

Given an elliptic curve point P of order n and a target point Q = kP, the algorithm generates a sequence of points Wi where each Wi can be expressed as Wi = aiP + biQ for known coefficients ai and bi.

The sequence is designed to eventually cycle, and when we detect Wi = Wj for i ≠ j, we have:
```
aiP + biQ = ajP + bjQ
(ai - aj)P = (bj - bi)Q
(ai - aj)P = (bj - bi)kP
```

Solving for k: `k ≡ (ai - aj)(bj - bi)^(-1) (mod n)`

### Walking Strategy

The algorithm partitions the elliptic curve group into r regions (default r=3) based on the x-coordinate:
- **Zone 0** (x ≡ 0 mod r): Wi+1 = Wi + P, ai+1 = ai + 1, bi+1 = bi
- **Zone 1** (x ≡ 1 mod r): Wi+1 = Wi + Q, ai+1 = ai, bi+1 = bi + 1  
- **Zone 2** (x ≡ 2 mod r): Wi+1 = 2Wi, ai+1 = 2ai, bi+1 = 2bi


## Complexity Analysis

### Time Complexity
- **Expected**: O(√n) where n is the order of point P
- **Theoretical bound**: ~1.25√(πn/2) group operations
- **Practical performance**: ??? TODO

### Space Complexity  
- **O(1)** - Only stores current state of tortoise and hare
- Additional O(k) for visualization traces (optional)


## Usage Examples

### Basic Usage
```python
from elliptic_curve import *
from basic_rho import *

# Define curve y² = x³ + 2x + 3 mod 97
curve = EllipticCurve(a=2, b=3, p=97)
P = Point(3, 6, curve)
order = curve.find_order(P)

# Create discrete logarithm instance
k_secret = 7
Q = curve.scalar_mul(k_secret, P)

# Solve using basic rho
k_found, _, _, steps = pollard_rho(P, Q, order, curve)
print(f"Found k = {k_found} in {steps} steps")
```

### With Visualization
```python
from basic_rho import *

# Run with step-by-step animation
k_found = animated_visualization(P, Q, order, curve, r=3)
```

## Visualization Features

The implementation includes an animated visualization that shows:
- **Tortoise path** (blue line): Single-step progression
- **Hare path** (red dashed line): Double-step progression  
- **Collision point** (gold star): Where tortoise and hare meet
- **Step counter**: Progress through the algorithm
- **Directional arrows**: Movement direction in each step