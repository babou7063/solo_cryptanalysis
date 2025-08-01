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

### Distinguished Points Strategy

The algorithm uses "distinguished points" : special points that satisfy a predetermined property (e.g., x-coordinate has specific bits). Walks stop when they reach a distinguished point, enabling:
- **Parallelization**: Multiple independent walks
- **Memory efficiency**: Store only distinguished points
- **Collision detection**: Compare distinguished points from different walks


## Complexity Analysis

### Time Complexity
- **Expected**: ??? TODO
- **Theoretical bound**: ??? TODO
- **Practical performance**: ??? TODO

### Space Complexity  
- ??? TODO


### Distinguished Point Density
The fraction of distinguished points affects performance:
- **Too dense** (e.g., 1/4): Short walks, high overhead
- **Too sparse** (e.g., 1/1024): Long walks, memory issues  
- **Optimal** (typically 1/16 to 1/64): Balance walk length and collision detection
