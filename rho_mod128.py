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
        current = self.P.copy()
        
        for i in range(self.r):
            # R_j = c_j * P + d_j * Q where c_j, d_j are random
            c_j = random.randint(1, self.order - 1)
            d_j = random.randint(1, self.order - 1)
            
            # Simple version with modular addition
            table.append((current.copy(), c_j, d_j))
            current = current.add(self.P)
        
        return table

    def hash_function(self, point: Point):
        
        # Use y-coordinate as int the paper
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
        if len(history) >= 4 and history[-1] == history[-3]:
            self.fruitless_cycles_detected += 1
            return True
        
        # Check for 4-cycle: W_{i-1} == W_{i-5}
        if len(history) >= 6 and history[-1] == history[-5]:
            self.fruitless_cycles_detected += 1
            return True
        
        return False

    def escape_fruitless_cycle(self, history: list):
        self.cycle_escapes += 1
        
        # Find minimum point in recent history
        recent_points = history[-4:] if len(history) >= 4 else history
        # Point with the min abs(x)
        min_point = min(recent_points, key=lambda p: p.x.to_int())
        
        # Double the point to get out of the cycle
        return min_point.double().canonical_form()

    def is_distinguished_point(self, point: Point) -> bool:
        # Simple criterion: x-coordinate has many trailing zeros
        x_int = point.x.to_int()
        return (x_int & 0xFFFFF) == 0  # 20 trailing zero bits
    
    def negative_additive_walk(self, max_iterations: int = 1000000):
        # Random starting point
        a = random.randint(1, self.order - 1)
        b = random.randint(1, self.order - 1)
        
        W = self.P.copy()  # Simplified starting point
        
        history = [W.copy()]
        cycle_check_frequency = max(1, int(2 * (self.r ** 0.5)))  # Check every ~2√r iterations
        
        for i in range(max_iterations):
            # Check for fruitless cycles periodically
            if i % cycle_check_frequency == 0 and self.detect_fruitless_cycle(history):
                W = self.escape_fruitless_cycle(history)
                # Restart with new starting point
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

    def solve_discrete_logarithm():
        #...
        return 9
    
    


# TEST
P = Point( Mod128Minus3Element.from_int(12345), Mod128Minus3Element.from_int(67890))
Q = Point( Mod128Minus3Element.from_int(54321), Mod128Minus3Element.from_int(98765))


ec = EllipticCurve(P, Q, 10)
print(ec.precomputed_table)