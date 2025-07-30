from mod128_minus3 import Mod128Minus3Element
import random

def extended_gcd(a, b):
    if a == 0:
        return b, 0, 1
    gcd, x1, y1 = extended_gcd(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return gcd, x, y

def modinv2(a, m):
    if a < 0:
        a = (a % m + m) % m
    g, x, _ = extended_gcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    return x % m

def simple_update(a, b, z, g, h, order, r=3):
    zone = z.to_int() % r
    
    if zone == 0:
        a_new = (a + 1) % order
        b_new = b
        z_new = (z * g).reduce()
    elif zone == 1:
        a_new = a
        b_new = (b + 1) % order
        z_new = (z * h).reduce()
    else:  # zone == 2
        a_new = (2 * a) % order
        b_new = (2 * b) % order
        z_new = z.square()
    return a_new, b_new, z_new

def pollard_rho_mod128_debug(g, h, order, r=3, max_iterations=10000):
    print(f"Starting Pollard rho with order={order}, max_iterations={max_iterations}")
    
    # Random starting point
    a = random.randint(0, order - 1)
    b = random.randint(0, order - 1)
    print(f"Starting with a={a}, b={b}")
    
    # Init z = g^a * h^b
    g_a = (g ** a).reduce()
    h_b = (h ** b).reduce()
    z = (g_a * h_b).reduce()
    print(f"Initial z = {z.to_int()}")

    # Init points of tortoise and hare
    A_tortoise, B_tortoise, Z_tortoise = a, b, z.copy()
    A_hare, B_hare, Z_hare = a, b, z.copy()

    seen_values = set()
    
    for i in range(max_iterations):
        # Tortoise: one step
        A_tortoise, B_tortoise, Z_tortoise = simple_update(A_tortoise, B_tortoise, Z_tortoise, g, h, order, r)
        
        # Hare: two steps
        A_hare, B_hare, Z_hare = simple_update(A_hare, B_hare, Z_hare, g, h, order, r)
        A_hare, B_hare, Z_hare = simple_update(A_hare, B_hare, Z_hare, g, h, order, r)
        
        # Debug: track values
        z_tort_int = Z_tortoise.to_int()
        z_hare_int = Z_hare.to_int()
        
        if i < 20:  # Print first 20 iterations
            print(f"Iter {i+1}: Tortoise z={z_tort_int}, Hare z={z_hare_int}")
        
        # Track seen values to detect if we're in a cycle
        if z_tort_int in seen_values:
            print(f"Tortoise revisited value {z_tort_int} at iteration {i+1}")
        seen_values.add(z_tort_int)
        
        # Check for collision
        if Z_tortoise == Z_hare:
            print(f"Collision found at iteration {i+1}")
            print(f"Tortoise: a={A_tortoise}, b={B_tortoise}, z={z_tort_int}")
            print(f"Hare: a={A_hare}, b={B_hare}, z={z_hare_int}")
            
            da = (A_tortoise - A_hare) % order
            db = (B_hare - B_tortoise) % order
            
            print(f"da={da}, db={db}")
            
            if db == 0:
                print("db = 0, continuing...")
                continue
            
            db_inv = modinv2(db, order)
            x = (da * db_inv) % order
            print(f"Computed x = {x}")
            return x
    
    print(f"No collision found. Saw {len(seen_values)} unique values.")
    print(f"Expected approximately sqrt(pi*{order}/2) = {(3.14159*order/2)**0.5:.1f} iterations")
    raise Exception("No collision found within iteration limit")

# Test 1: Very basic verification
print("=== Basic Verification ===")
q = 13
g = Mod128Minus3Element.from_int(2)
x_secret = 5
h = (g ** x_secret).reduce()

print(f"g = {g.to_int()}")
print(f"h = {h.to_int()}")
print(f"g^{x_secret} mod {q} = {pow(2, x_secret, q)}")
print(f"Verification: 2^5 = 32, 32 mod 13 = {32 % 13}")

# Test 2: Try with even smaller numbers
print("\n=== Smaller Test ===")
q_small = 7
g_small = Mod128Minus3Element.from_int(3)
x_small = 2
h_small = (g_small ** x_small).reduce()

print(f"Small test: g={g_small.to_int()}, x={x_small}, h={h_small.to_int()}")
print(f"Verification: 3^2 = 9, 9 mod 7 = {9 % 7}")

try:
    x_found = pollard_rho_mod128_debug(g_small, h_small, q_small, r=3, max_iterations=1000)
    print(f"SUCCESS! Found x = {x_found}")
except Exception as e:
    print(f"Failed: {e}")