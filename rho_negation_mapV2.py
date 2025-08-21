# Implementation adapted from:
# D. J. Bernstein, T. Lange, P. Schwabe, 
# "On the correct use of the negation map in the Pollard rho method", 2010.
# https://dl.acm.org/doi/10.5555/1964658.1964669
#
# This code follows the techniques described in the paper, particularly:
# - use of the negation map to reduce expected runtime by √2
# - additive walks with precomputed tables R_j = c_j*P + d_j*Q
# - detection and escape from fruitless cycles (length 2, 4, etc.)
# - branchless iteration and use of canonical point forms.
# 
# This adaptation is for educational and experimental purposes.

from elliptic_curve import Point, EllipticCurve
from utils import modinv
import random
import hashlib
import time
import math

class NegationMapRho:
    """
    Pollard rho for elliptic curves over small prime fields.
    Uses the negation map and additive walks like the mod128 version but simpler.
    """
    
    def __init__(self, curve, P, Q, order, r=64, dp_bits=None, alpha=8, max_iter_factor=8):
        """
        Initialize a NegationMapRho instance.

        :param curve: the elliptic curve
        :param P: the point P
        :param Q: the point Q
        :param order: the order of the curve
        :param r: (optional) the size of the precomputed table (default 64)
        """
        self.curve = curve
        self.P = P
        self.Q = Q
        self.order = order
        self.r = r
        
        self.curve.reset_op_counters()
        self.precomputed_table = self.generate_precomputed_table()
        self.precomp_ops = self.curve.get_op_counters()

        # --- dynamic DP ---
        if dp_bits is None:
            t = int(round(math.log2(max(2, int(math.isqrt(self.order) // max(1, alpha))))))
            self.dp_bits = max(1, min(24, t))
        else:
            self.dp_bits = max(1, int(dp_bits))
        self.dp_mask = (1 << self.dp_bits) - 1

        self.max_iterations = max(512, int(max_iter_factor * math.isqrt(self.order)))

        # Statistics
        self.fruitless_cycles_detected = 0
        self.cycle_escapes = 0
        self.walks_performed = 0
        
    def generate_precomputed_table(self):
        """
        Generate precomputed table R_j = c_j * P + d_j * Q.
        
        :return: a list of r tuples (Rj, cj, dj)
        """
        table = []
        
        i = 0
        while i < self.r:
            # Generate random coefficients
            c_j = random.randint(1, self.order - 1)
            d_j = random.randint(1, self.order - 1)
            
            # R_j = c_j * P + d_j * Q
            R_j = self.curve.scalar_mul(c_j, self.P) + self.curve.scalar_mul(d_j, self.Q)
            
            # Avoid point at infinity in table, if point at infinity just try again with different coefficients
            if not R_j.at_infinity:
                table.append((self.canonical_form(R_j), c_j, d_j))
                i += 1
        return table
    
    def canonical_form(self, point):
        """
        Return the canonical form of a point on the elliptic curve.
        The function selects between the point and its negation based on the
        parity of the y-coordinate.

        :param point: The point on the elliptic curve to convert to canonical form.
        :return: The canonical form of the point.
        """
        if point.at_infinity:
            return point
        
        # Choose point with even y-coordinate
        if point.y % 2 == 0:
            return point
        else:
            return Point(point.x, (-point.y) % self.curve.p, self.curve)
        
    def hash_function(self, point):
        """
        Hash function for table lookup.
        The function takes a point and returns an integer in the range [0, r-1].
        
        :param point: The point to hash.
        :return: An integer hash value.
        """
        if point.at_infinity:
            return 0
        
        # Use canonical form and hash y-coordinate
        canonical = self.canonical_form(point)
        return canonical.y % self.r
    
    def additive_walk_step(self, W, a, b):
        """
        Perform a single step of the additive walk in the Pollard's rho algorithm.
        This function updates the current point W on the elliptic curve and its
        associated coefficients (a, b) using a precomputed table of random linear
        combinations.

        :param W: The current point on the elliptic curve.
        :param a: The current coefficient of the base point P.
        :param b: The current coefficient of the target point Q.
        :return: A tuple containing the new canonical point on the curve and the
                updated coefficients (new_W_canonical, new_a, new_b).
        """

        canonical_W = self.canonical_form(W)
        hash_val = self.hash_function(canonical_W)
        
        R_h, c_h, d_h = self.precomputed_table[hash_val]
        
        # W_{i+1} = W_i + R_h
        new_W = canonical_W + R_h
        new_W_canonical = self.canonical_form(new_W)
        
        # Update coefficients
        new_a = (a + c_h) % self.order
        new_b = (b + d_h) % self.order
        
        return new_W_canonical, new_a, new_b
    
    def points_equal(self, P1, P2):
        """
        Check if two points on an elliptic curve are equal.

        :param P1: First point on the elliptic curve.
        :param P2: Second point on the elliptic curve.
        :return: True if the points are equal, False otherwise.
        """

        if P1.at_infinity and P2.at_infinity:
            return True
        if P1.at_infinity or P2.at_infinity:
            return False
        return P1.x == P2.x and P1.y == P2.y

    def detect_fruitless_cycle(self, history):
        """
        Detect fruitless cycles in the Pollard's rho algorithm.

        :param history: A list of points visited during the walk.
        :return: True if a fruitless cycle has been detected, False otherwise.
        """
        
        if len(history) < 4:
            return False
        
        # Check for 2-cycle: W_{i-1} == W_{i-3}
        if len(history) >= 4:
            if self.points_equal(history[-1], history[-3]):
                self.fruitless_cycles_detected += 1
                return True
        
        # Check for 4-cycle: W_{i-1} == W_{i-5}
        if len(history) >= 6:
            if self.points_equal(history[-1], history[-5]):
                self.fruitless_cycles_detected += 1
                return True
        
        return False
    
    
    def escape_fruitless_cycle(self, history):
        """
        Escape from a detected fruitless cycle by doubling the minimum point.

        :param history: A list of points visited during the walk.
        :return: A point that is the canonical form of the doubled minimum point.
        """

        self.cycle_escapes += 1
        
        # Find minimum point in recent history (by x-coordinate)
        recent_points = history[-4:] if len(history) >= 4 else history
        min_point = min(recent_points, key=lambda p: p.x if not p.at_infinity else float('inf'))
        
        # Double the point to escape
        return self.canonical_form(min_point + min_point)
    
    def is_distinguished_point(self, point):
        """
        Determine if a given point is a distinguished point in the context of the Pollard's rho algorithm.
        A point is considered distinguished if it is not at infinity and its
        x-coordinate has a specific number of trailing zero bits.

        :param point: The point on the elliptic curve to check.
        :return: True if the point is distinguished, False otherwise.
        """

        if point.at_infinity:
            return False
        return (point.x & self.dp_mask) == 0
    
    def single_walk(self, max_iterations=None):
        """
        Perform a single walk until hitting a distinguished point.
        Start at a random point W = aP + bQ
        Then add a point Rj from a precomputed table to W
        Stops when it finds a distinguished point

        :param max_iterations: The maximum number of steps to take in the walk.
        :return: A tuple of the distinguished point, the final coefficients (a, b), and the number of steps taken.
        """
        
        # Random starting coefficients
        a = random.randint(1, self.order - 1)
        b = random.randint(1, self.order - 1)
        
        # Starting point W = a*P + b*Q
        W = self.curve.scalar_mul(a, self.P) + self.curve.scalar_mul(b, self.Q)
        W = self.canonical_form(W)
        
        history = [W]
        cycle_check_freq = max(1, int(2 * (self.r ** 0.5)))  # Check every ~2√r steps

        limit = self.max_iterations if max_iterations is None else max_iterations
        for i in range(limit):
            # Periodically check for fruitless cycles
            if i % cycle_check_freq == 0 and self.detect_fruitless_cycle(history):
                W = self.escape_fruitless_cycle(history)
                # Reset coefficients after escape
                a = random.randint(1, self.order - 1)
                b = random.randint(1, self.order - 1)
                history = [W]
                continue
            
            # Perform additive walk step
            W, a, b = self.additive_walk_step(W, a, b)
            history.append(W)
            
            # Keep history bounded
            if len(history) > 10:
                history.pop(0)
            
            # Check for distinguished point
            if self.is_distinguished_point(W):
                return W, a, b, i + 1
        
        return None
    
    def solve_ecdlp(self, max_walks=1000, return_stats=False):
        """
        Solve the elliptic curve discrete logarithm problem using the Pollard rho algorithm.

        This method starts at a random point W = aP + bQ on the elliptic curve and iteratively
        adds a point Rj from a precomputed table to W until it finds a distinguished point.
        The algorithm then checks for collisions against previously seen distinguished points
        and if a collision is found, solves for the discrete logarithm.

        :param max_walks: The maximum number of walks to perform before giving up.
        :return: The discrete logarithm k such that Q = kP, or None if no solution is found
        """
        #print(f"Starting small-scale Pollard rho with:")
        #print(f"  Curve: y² = x³ + {self.curve.a}x + {self.curve.b} mod {self.curve.p}")
        #print(f"  Order: {self.order}")
        #print(f"  Table size r: {self.r}")
        #print(f"  Expected ~{(3.14159 * self.order / 4) ** 0.5:.0f} iterations per walk")
        
        self.curve.reset_op_counters()
        distinguished_points = {}
        
        for walk_id in range(max_walks):
            self.walks_performed += 1
            
            #if walk_id % 50 == 0 and walk_id > 0:
                #print(f"Walk {walk_id}, cycles: {self.fruitless_cycles_detected}, escapes: {self.cycle_escapes}")
            
            # Perform single walk
            result = self.single_walk()
            if result is None:
                continue
            
            W, a, b, iterations = result
            
            # Create hashable key for the point
            point_key = (W.x, W.y) if not W.at_infinity else ("inf", "inf")
            
            # Check for collision
            if point_key in distinguished_points:
                a_prev, b_prev, prev_walk = distinguished_points[point_key]
                
                #print(f"\nCollision found between walks {prev_walk} and {walk_id}!")
                #print(f"Distinguished point: {W}")
                #print(f"Walk lengths: {iterations} iterations")
                
                # Solve for discrete logarithm
                da = (a - a_prev) % self.order
                db = (b_prev - b) % self.order
                
                if db == 0:
                    #print("db = 0, continuing...")
                    continue
                
                try:
                    db_inv = modinv(db, self.order)
                    k = (da * db_inv) % self.order
                    
                    #print(f"Found discrete logarithm: k = {k}")
                    
                    # Verify the result
                    verification = self.curve.scalar_mul(k, self.P)
                    if verification == self.Q:
                        walk_ops = self.curve.get_op_counters()
                        stats = {
                            "ops_adds": walk_ops["adds"] + self.precomp_ops["adds"],
                            "ops_dbls": walk_ops["dbls"] + self.precomp_ops["dbls"],
                            "ops_total": walk_ops["total"] + self.precomp_ops["total"],
                            "ops_precompute": self.precomp_ops["total"],
                            "walks": self.walks_performed,
                            "r": self.r
                        }
                        return k, stats
                    else:
                        continue
                        
                except Exception as e:
                    continue
            else:
                distinguished_points[point_key] = (a, b, walk_id)
        return None


# TEST

def run_simple_tests():
    from time import time

    for p, a, b, k, label in [(101, 2, 3, 7, "p=101"), (1009, 2, 3, 123, "p=1009")]:
        curve = EllipticCurve(a, b, p)
        P = next(Point(x, y, curve) for x in range(p) for y in range(p)
                 if (y*y) % p == (x**3 + a*x + b) % p)
        Q = curve.scalar_mul(k, P)
        order = curve.find_order(P)

        print(f"\n=== Test: {label} ===")
        print(f"P = {P}, Q = {Q}, k = {k}, order = {order}")

        solver = NegationMapRho(curve, P, Q, order, r=32)
        t0 = time()
        found_k = solver.solve_ecdlp(max_walks=500)
        t1 = time()

        if found_k is not None and (found_k - k) % order == 0:
            print(f"Found k = {found_k} in {t1 - t0:.2f}s")
        else:
            print(f"Failed (k = {found_k}) in {t1 - t0:.2f}s")

#run_simple_tests()