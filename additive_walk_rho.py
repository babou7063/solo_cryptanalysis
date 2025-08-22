from elliptic_curve import Point, EllipticCurve
from utils import modinv
import random
import hashlib
import math

def is_distinguished(W):
    """
    Determine if a given point W is a "distinguished" point.
    A point is considered distinguished if it is not the point at infinity
    and its x-coordinate is congruent to 0 modulo 32.

    :param W: The point on the elliptic curve to check
    :return: True if the point is distinguished, False otherwise
    """
    return not W.at_infinity and W.x % 8 == 0


def precompute_table(P, Q, r, order, curve):
    """
    Precompute a lookup table for the additive walk rho algorithm.

    Given points P and Q on a curve, and parameters r and order, this function
    precomputes a lookup table of r points Rj, along with coefficients cj and dj
    such that Rj = cjP + djQ, for 0 <= j < r. The coefficients cj and dj are
    randomly chosen from the range [0, order - 1].

    :param P: a point on the curve
    :param Q: a point on the curve
    :param r: the number of points in the lookup table
    :param order: the range of the coefficient values
    :param curve: the elliptic curve
    :return: a list of r tuples (Rj, cj, dj)
    """
    
    table = []
    for _ in range(r):
        while True:
            cj = random.randint(0, order - 1)
            dj = random.randint(0, order - 1)
            if cj == 0 and dj == 0:
                continue
            Rj = curve.scalar_mul(cj, P) + curve.scalar_mul(dj, Q)
            if not Rj.at_infinity:
                break
        table.append((Rj, cj, dj))
    return table


def h(point, r):
    """
    Compute a hash value for a given point on the elliptic curve.
    This function first converts the x and y coordinates of the point to bytes
    and then computes the SHA-256 hash of those bytes. The first byte of the
    resulting hash is then returned modulo r.

    :param point: The point on the elliptic curve to hash
    :param r: The modulus to reduce the hash value to
    :return: The hash value of the point modulo r
    """
    if point.at_infinity:
        return 0
    x = point.x
    y = point.y
    xb = x.to_bytes((x.bit_length() + 7)//8 or 1, 'big')
    yb = y.to_bytes((y.bit_length() + 7)//8 or 1, 'big')
    digest = hashlib.sha256(b'X'+xb+b'Y'+yb).digest()
    
    return digest[0] % r


def additive_walk(P, Q, order, curve, r=16, max_steps=10_000, is_distinguished=None, R_table=None):
    """
    Run the additive walk rho algorithm to find a collision between a "distinguished"
    point and a random walk on the elliptic curve.

    The algorithm starts at a random point W = aP + bQ on the curve, and then
    iteratively adds a point Rj from a precomputed table to W. The index j is
    chosen by hashing the current value of W modulo r.

    The algorithm stops when it finds a "distinguished" point.
    The predicate is passed as the is_distinguished
    argument, and defaults to a simple predicate that holds when the x-coordinate
    of the point is congruent to 0 modulo 32.

    :param P: a point on the elliptic curve
    :param Q: a point on the elliptic curve
    :param order: the range of the coefficient values
    :param curve: the elliptic curve
    :param r: the number of points in the lookup table
    :param max_steps: the maximum number of steps to take in the walk
    :param is_distinguished: a predicate that takes a point on the curve as input
        and returns a boolean indicating whether the point is distinguished
    :return: a tuple of the distinguished point, the random seed used to start
        the walk, and the precomputed table of points
    """
    if is_distinguished is None:
        is_distinguished = lambda W: (not W.at_infinity) and (W.x % 32 == 0)

    if R_table is None:
        R_table = precompute_table(P, Q, r, order, curve)

    # Random seed: a, b
    a0 = random.randint(1, order - 1)
    b0 = random.randint(1, order - 1)
    seed = (a0, b0)

    W = curve.scalar_mul(a0, P) + curve.scalar_mul(b0, Q)

    for step in range(max_steps):
        if is_distinguished(W):
            return W, seed, R_table

        j = h(W, r)
        Rj, cj, dj = R_table[j]
        W = W + Rj

    raise Exception("No distinguished point found")


def replay_walk(seed, R_table, P, Q, order, curve, target_point, r, max_steps=10_000):
    
    """
    Replay the additive walk to find coefficients a and b such that target_point = aP + bQ.
    Stop when W == target_point and return a, b.
    
    :param seed: a tuple (a, b) representing the initial coefficients
    :param R_table: a precomputed lookup table of points and coefficients
    :param P: a point on the elliptic curve
    :param Q: a point on the elliptic curve
    :param order: the range of the coefficient values
    :param curve: the elliptic curve
    :param target_point: the point to reach during the walk
    :param r: the number of points in the lookup table
    :return: a tuple (a, b) of the coefficients such that target_point = aP + bQ
    """

    a, b = seed
    W = curve.scalar_mul(a, P) + curve.scalar_mul(b, Q)

    for step in range(max_steps):
        if W == target_point:
            return a, b
        j = h(W, r)
        Rj, cj, dj = R_table[j]
        W = W + Rj
        a = (a + cj) % order
        b = (b + dj) % order

    raise Exception("Replay walk exceeded max_steps - possible desynchronization")

def solve_for_k(num, den, n, P, Q, curve):
    
    """
    Solve for k in the equation num*k = den (mod n), where n is the order of P.
    If there is no solution, return None.

    :param num: The left-hand side of the equation
    :param den: The right-hand side of the equation
    :param n: The order of P
    :param P: The point on the elliptic curve
    :param Q: The point on the elliptic curve
    :param curve: The elliptic curve
    :return: The solution k if it exists, None otherwise
    """
    d = math.gcd(den, n)
    if num % d != 0:
        return None
    n1, den1, num1 = n // d, den // d, num // d
    try:
        # Compute the modular inverse of den1 modulo n1
        inv = pow(den1, -1, n1) 
    except ValueError:
        return None
    k0 = (num1 * inv) % n1
    for t in range(d):
        k = (k0 + t * n1) % n
        if curve.scalar_mul(k, P) == Q:
            return k
    return None

def retry_walks(P, Q, order, curve, r, is_distinguished, max_attempts=10, walks_per_attempt=64, max_steps=10_000, return_stats=False):
    
    for attempt in range(max_attempts):
        curve.reset_op_counters()
        R_table = precompute_table(P, Q, r, order, curve)
        ops_after_precomp = curve.get_op_counters()
        precomp_ops = ops_after_precomp["total"]

        seen = {}
        walks = 0
        replays = 0

        for _ in range(walks_per_attempt):
            walks += 1
            try:
                W, seed, _ = additive_walk(P, Q, order, curve, r=r, max_steps=max_steps, is_distinguished=is_distinguished, R_table=R_table)
            except Exception:
                continue

            key = (W.x, W.y) if not W.at_infinity else ('INF', 'INF')

            if key in seen:
                seed1 = seen[key]
                seed2 = seed
                replays += 2
                try:
                    a1, b1 = replay_walk(seed1, R_table, P, Q, order, curve, W, r, max_steps=max_steps)
                    a2, b2 = replay_walk(seed2, R_table, P, Q, order, curve, W, r, max_steps=max_steps)
                except Exception:
                    continue

                num = (a1 - a2) % order
                den = (b2 - b1) % order
                k = solve_for_k(num, den, order, P, Q, curve)
                if k is None:
                    continue

                ops_now = curve.get_op_counters()
                stats = {
                    "ops_adds": ops_now["adds"],
                    "ops_dbls": ops_now["dbls"],
                    "ops_total": ops_now["total"],
                    "ops_precompute": precomp_ops,
                    "walks": walks,
                    "replays": replays,
                    "r": r
                }
                return (k, stats) if return_stats else k
            else:
                seen[key] = seed

    return (None, None) if return_stats else None

"""
##########################################
#                  Test                  #
##########################################
p = 211
a = 2
b = 3
curve = EllipticCurve(a, b, p)

P = Point(1, 39, curve)
order = curve.find_order(P)
print(f"Order of P: {order}")

k_secret = 7
Q = curve.scalar_mul(k_secret, P)

k_found = retry_walks(P, Q, order, curve, r=8, is_distinguished=is_distinguished, max_attempts=10, walks_per_attempt=64, max_steps=10_000)
print(f"Recovered k = {k_found}")
"""
