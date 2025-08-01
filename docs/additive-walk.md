# Additive Walk Pollard's Rho Algorithm

## Overview

The Additive Walk variant of Pollard's rho algorithm improves upon the basic version by using precomputed tables to create more random walks and enable better parallelization through distinguished points. This approach typically achieves better performance and is more suitable for distributed implementations.

## Algorithm Description

### Core Innovation

Instead of using simple partition-based rules, the additive walk uses a precomputed table of random linear combinations. Each step adds a table entry based on a hash of the current point:

```
Wi+1 = Wi + Rj    where j = h(Wi) mod r
Rj = cj·P + dj·Q  (precomputed with random cj, dj)
```

The coefficients are tracked as:
```
ai+1 = ai + cj (mod order)
bi+1 = bi + dj (mod order)
```
