from elliptic_curve import Point, EllipticCurve
from utils import modinv
import random
import hashlib
import time

class NegationMapRho:
    """
    Pollard rho for elliptic curves over small prime fields.
    Uses the negation map and additive walks like the mod128 version but simpler.
    """
    
    def __init__(self, curve, P, Q, order, r=64):
        """
        Initialize the small-scale ECDLP solver.
        
        :param curve: EllipticCurve instance
        :param P: Base point
        :param Q: Target point (Q = k*P, we want to find k)
        :param order: Order of point P
        :param r: Size of precomputed table
        """
        self.curve = curve
        self.P = P
        self.Q = Q
        self.order = order
        self.r = r
        self.precomputed_table = self.generate_precomputed_table()
        
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
    
    def detect_fruitless_cycle(self, history):
        return 0
    
    
    def escape_fruitless_cycle(self, history):
        return 0
    
    def is_distinguished_point(self, point):
        return 0
    
    def single_walk(self, max_iterations=100000):
        return 0
    
    def solve_ecdlp(self, max_walks=1000):
        return 0

