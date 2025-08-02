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

