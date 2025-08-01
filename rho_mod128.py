# Implementation adapted from:
# D. J. Bernstein, T. Lange, P. Schwabe, 
# "On the correct use of the negation map in the Pollard rho method", 2010.
# https://dl.acm.org/doi/10.5555/1964658.1964669
#
# This code implements key ideas from the paper, including:
# - the use of the negation map to reduce the expected runtime of Pollard rho
# - additive walks with precomputed tables R_j = c_j·P + d_j·Q
# - canonical forms of points to avoid duplicates due to symmetry
# - detection and escape from fruitless cycles
# - optimized modular arithmetic for 2^128 - 3
# 
# This adaptation is for educational and experimental purposes.


from mod128_minus3 import Mod128Minus3Element, Point
import random

class EllipticCurve:
    
    def __init__(self, base_point: Point, target_point: Point, order: int, r: int = 2048):
        self.P = base_point  # Base point
        self.Q = target_point  # Target point (Q = kP, we want to find k)
        self.order = order
        self.r = r  # Table size
        self.precomputed_table = self.generate_precomputed_table()
        
        # For statistics
        self.fruitless_cycles_detected = 0
        self.cycle_escapes = 0
        
    def generate_precomputed_table(self):
        table = []
        
        for i in range(self.r):
            # R_j = c_j * P + d_j * Q where c_j, d_j are random
            c_j = random.randint(1, self.order - 1)
            d_j = random.randint(1, self.order - 1)
            
            # Compute R_j = c_j * P + d_j * Q
            R_j = self.scalar_mult(self.P, c_j).add(self.scalar_mult(self.Q, d_j))
            
            table.append((R_j.canonical_form(), c_j, d_j))
        
        return table

    def scalar_mult(self, point: Point, scalar: int) -> Point:
        """double-and-add"""
        if scalar == 0:
            result = Point.at_infinity()  # Point at infinity
        if scalar == 1:
            return point.copy()
        
        result = Point.at_infinity()  # Point at infinity
        addend = point.copy()
        
        while scalar > 0:
            if scalar & 1:
                if self.is_point_at_infinity(result):
                    result = addend.copy()
                else:
                    result = result.add(addend)
            addend = addend.double()
            scalar >>= 1
        
        return result

    def is_point_at_infinity(self, point: Point) -> bool:
        """Check if point is at infinity (represented as (0,0))"""
        return point.x.to_int() == 0 and point.y.to_int() == 0

    def hash_function(self, point: Point):
        # Use y-coordinate as in the paper
        canonic = point.canonical_form()
        y_int = canonic.y.to_int()
        return y_int % self.r

    def additive_walk(self, W: Point, a: int, b: int):
        canonic_W = W.canonical_form()
        hash_W = self.hash_function(canonic_W)
        
        R_h, c_h, d_h = self.precomputed_table[hash_W]
        
        # W_{i+1} = W_i + R_h
        new_W = canonic_W.add(R_h)
        new_W_canonic = new_W.canonical_form()
        
        # Update coefficients
        new_a = (a + c_h) % self.order
        new_b = (b + d_h) % self.order
        
        return new_W_canonic, new_a, new_b

    def detect_fruitless_cycle(self, history: list):
        # Not enough history
        if len(history) < 4:
            return False
        
        # Check for 2-cycle: W_{i-1} == W_{i-3}
        if len(history) >= 4:
            p1, p3 = history[-1], history[-3]
            if p1.x.to_int() == p3.x.to_int() and p1.y.to_int() == p3.y.to_int():
                self.fruitless_cycles_detected += 1
                return True
        
        # Check for 4-cycle: W_{i-1} == W_{i-5}
        if len(history) >= 6:
            p1, p5 = history[-1], history[-5]
            if p1.x.to_int() == p5.x.to_int() and p1.y.to_int() == p5.y.to_int():
                self.fruitless_cycles_detected += 1
                return True
        
        return False

    def escape_fruitless_cycle(self, history: list):
        self.cycle_escapes += 1
        
        # Find minimum point in recent history
        recent_points = history[-4:] if len(history) >= 4 else history
        # Point with the min abs(x)
        min_point = min(recent_points, key=lambda p: abs(p.x.to_int()))
        
        # Double the point to get out of the cycle
        return min_point.double().canonical_form()

    def is_distinguished_point(self, point: Point) -> bool:
        # Simple criterion: x-coordinate has many trailing zeros
        # Skip point at infinity
        if self.is_point_at_infinity(point):
            return False
        x_int = point.x.to_int()
        return (x_int & 0xF) == 0  # 4 trailing zero bits (reduced for faster convergence)
    
    def negative_additive_walk(self, max_iterations: int = 1000000):
        # Random starting point: W = a*P + b*Q
        a = random.randint(1, self.order - 1)
        b = random.randint(1, self.order - 1)
        
        W = self.scalar_mult(self.P, a).add(self.scalar_mult(self.Q, b))
        W = W.canonical_form()
        
        history = [W.copy()]
        cycle_check_frequency = max(1, int(2 * (self.r ** 0.5)))  # Check every ~2√r iterations
        
        for i in range(max_iterations):
            # Check for fruitless cycles periodically
            if i % cycle_check_frequency == 0 and self.detect_fruitless_cycle(history):
                W = self.escape_fruitless_cycle(history)
                # Update coefficients randomly after escape
                a = random.randint(1, self.order - 1)
                b = random.randint(1, self.order - 1)
                history = [W.copy()]
                continue
            
            # Perform additive step
            W, a, b = self.additive_walk(W, a, b)
            history.append(W.copy())
            
            # Keep history bounded
            if len(history) > 10:
                history.pop(0)
            
            # Check if this is a distinguished point
            if self.is_distinguished_point(W):
                return W, a, b, i + 1
        
        return None

    def solve_discrete_logarithm(self, max_walks: int = 10000):
        print(f"Starting negating Pollard rho with r={self.r}")
        print(f"Expected iterations per walk: ~{(3.14159 * self.order / 4) ** 0.5:.0f}")
        
        distinguished_points = {}
        
        for walk_id in range(max_walks):
            if walk_id % 100 == 0:
                print(f"Walk {walk_id}, cycles detected: {self.fruitless_cycles_detected}, escapes: {self.cycle_escapes}")
            
            # Perform negative additive walk
            result = self.negative_additive_walk()
            if result is None:
                continue
            
            W, a, b, iterations = result
            
            # Convert point to a hashable representation
            point_key = (W.x.to_int(), W.y.to_int())
            
            # Check for collision
            if point_key in distinguished_points:
                a_prev, b_prev, prev_walk = distinguished_points[point_key]
                
                print(f"Collision found between walks {prev_walk} and {walk_id}")
                print(f"Point: x={W.x.to_int()}, y={W.y.to_int()}")
                
                # Solve for discrete logarithm
                da = (a - a_prev) % self.order
                db = (b_prev - b) % self.order
                
                if db == 0:
                    print("db = 0, continuing...")
                    continue
                
                try:
                    # Use extended Euclidean algorithm for modular inverse
                    def mod_inverse(a, m):
                        if m == 1:
                            return 0
                        m0, x0, x1 = m, 0, 1
                        while a > 1:
                            q = a // m
                            m, a = a % m, m
                            x0, x1 = x1 - q * x0, x0
                        return x1 % m0
                    
                    db_inv = mod_inverse(db, self.order)
                    k = (da * db_inv) % self.order
                    
                    print(f"Found discrete logarithm: k = {k}")
                    print(f"Verification: computing k*P...")
                    
                    # Verify the result
                    verification = self.scalar_mult(self.P, k)
                    if (verification.x.to_int() == self.Q.x.to_int() and 
                        verification.y.to_int() == self.Q.y.to_int()):
                        print("Verification successful!")
                    else:
                        print("Verification failed!")
                    
                    print(f"Total fruitless cycles detected: {self.fruitless_cycles_detected}")
                    print(f"Total cycle escapes: {self.cycle_escapes}")
                    
                    return k
                except Exception as e:
                    print(f"Failed to compute solution: {e}")
                    continue
            else:
                distinguished_points[point_key] = (a, b, walk_id)
        
        print("No collision found within walk limit")
        return None


# CORRECTED TEST with proper relationship Q = k*P
def create_test_case():
    # Create a base point P
    P = Point(Mod128Minus3Element.from_int(12345), Mod128Minus3Element.from_int(67890))
    
    # Choose a secret k
    k_secret = 157
    order = 1009
    
    # Create solver instance to use scalar_mult
    temp_solver = EllipticCurve(P, P, order, r=32)  # temporary
    
    # Compute Q = k*P
    Q = temp_solver.scalar_mult(P, k_secret)
    
    print(f"Created test case:")
    print(f"P = ({P.x.to_int()}, {P.y.to_int()})")
    print(f"Q = ({Q.x.to_int()}, {Q.y.to_int()})")
    print(f"k_secret = {k_secret}")
    
    return P, Q, k_secret

# TEST




"""
# Create a base point P
P = Point(Mod128Minus3Element.from_int(12345), Mod128Minus3Element.from_int(67890))
print(f"P = ({P.x.to_int()}, {P.y.to_int()})")

P2 = P.double()
print(f"P2 = ({P2.x.to_int()}, {P2.y.to_int()})")

P_add = P.add(P)
print(f"P + P = ({P_add.x.to_int()}, {P_add.y.to_int()})")

if P2 == P_add:
    print("\n TEST PASSED: P.double() == P.add(P)")
else:
    print("\n TEST FAILED: P.double() != P.add(P)")
"""

P, Q, k_secret = create_test_case()
order = 1009

solver = EllipticCurve(P, Q, order, r=32)

print(f"\nSearching for k such that Q = k*P with order {order}")
print(f"Expected k = {k_secret}")

found_k = solver.solve_discrete_logarithm(max_walks=1000)

if found_k is not None:
    print(f"SUCCESS: Found k = {found_k}")
    if (found_k - k_secret) % order == 0:
        print("SUCCESS: Correct discrete logarithm found!")
    else:
        print("FAILED: k is not correct modulo order")
        print(f"Difference: {(found_k - k_secret) % order}")
else:
    print("FAILED: No discrete logarithm found")
    