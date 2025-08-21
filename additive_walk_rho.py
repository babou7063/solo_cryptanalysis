from elliptic_curve import Point, EllipticCurve
from utils import modinv
import random
import hashlib

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


def additive_walk(P, Q, order, curve, r=16, max_steps=10_000, is_distinguished=None):
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
        is_distinguished = lambda W: W.x % 32 == 0  # Ex: x ≡ 0 mod 32

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
        W = W + R_table[j][0]

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


def retry_walks(P, Q, order, curve, r, is_distinguished, max_attempts=10):
    """
    Attempt to recover the scalar k such that Q = kP using the additive walk rho algorithm.

    :param P: a point on the elliptic curve
    :param Q: a point on the elliptic curve, such that Q = kP
    :param order: the order of the group generated by P
    :param curve: the elliptic curve
    :param r: the number of partitions in the walk
    :param is_distinguished: a predicate function to determine if a point is "distinguished"
    :param max_attempts: maximum number of attempts to find a collision (default: 10)
    :return: the recovered scalar k, or None if all attempts fail
    """

    for attempt in range(max_attempts):
        try:
            W1, seed1, R_table = additive_walk(P, Q, order, curve, r=r, is_distinguished=is_distinguished)
            W2, seed2, _ = additive_walk(P, Q, order, curve, r=r, is_distinguished=is_distinguished)

            if W1 == W2 and seed1 == seed2:
                continue  # same seed ⇒ unuseful

            if W1 == W2:
                #print("Collision detected!")
                a1, b1 = replay_walk(seed1, R_table, P, Q, order, curve, W1, r)
                a2, b2 = replay_walk(seed2, R_table, P, Q, order, curve, W2, r)

                if (b2 - b1) % order == 0:
                    #print("Failure: division by 0")
                    continue
                k = ((a1 - a2) * modinv(b2 - b1, order)) % order
                return k
        except Exception as e:
            print(f"Attempt {attempt + 1} failed: {e}")
    
    #print("All attempts failed.")
    return None


##########################################
#                  Test                  #
##########################################
"""
p = 211
a = 2
b = 3
curve = EllipticCurve(a, b, p)

P = Point(1, 39, curve)
order = curve.find_order(P)
print(f"Order of P: {order}")

k_secret = 7
Q = curve.scalar_mul(k_secret, P)

k_found = retry_walks(P, Q, order, curve, r=8, is_distinguished=is_distinguished)
print(f"Recovered k = {k_found}")


"""