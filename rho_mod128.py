from mod128_minus3 import Mod128Minus3Element
import random
from utils import modinv2

def partition(z, r=7):
    return z.to_int() % r

def update(a, b, z, g, h, order, r=7):
    zone = partition(z, r)
    if zone == 0:
        # z := z * g, a := a + 1, b unchanged
        a_new = (a + 1) % order
        b_new = b
        z_new = (z * g).reduce()
    elif zone == 1:
        # z := z * h, b := b + 1, a unchanged
        a_new = a
        b_new = (b + 1) % order
        z_new = (z * h).reduce()
    elif zone == 2:
        # z := z^2, a := 2*a, b := 2*b
        a_new = (2 * a) % order
        b_new = (2 * b) % order
        z_new = z.square().reduce()
    elif zone == 3:
        # z := z * g^2, a := a + 2, b unchanged
        a_new = (a + 2) % order
        b_new = b
        z_new = (z * g * g).reduce()
    elif zone == 4:
        # z := z * h^2, b := b + 2, a unchanged
        a_new = a
        b_new = (b + 2) % order
        z_new = (z * h * h).reduce()
    elif zone == 5:
        # z := z^3, a := 3*a, b := 3*b
        a_new = (3 * a) % order
        b_new = (3 * b) % order
        z_new = (z * z * z).reduce()
    else:
        # z := z * g * h, a := a + 1, b := b + 1
        a_new = (a + 1) % order
        b_new = (b + 1) % order
        gh = (g * h).reduce()
        z_new = (z * gh).reduce()

    return a_new, b_new, z_new

def simple_update(a, b, z, g, h, order, r=3):
    """
    Simplified update function similar to basic elliptic curve rho.
    This should provide better mixing than the complex 7-way partition.
    """
    zone = partition(z, r)
    
    if zone == 0:
        # z := z * g, a := a + 1
        a_new = (a + 1) % order
        b_new = b
        z_new = (z * g).reduce()
    elif zone == 1:
        # z := z * h, b := b + 1
        a_new = a
        b_new = (b + 1) % order
        z_new = (z * h).reduce()
    else:  # zone == 2
        # z := z^2, a := 2*a, b := 2*b
        a_new = (2 * a) % order
        b_new = (2 * b) % order
        z_new = z.square()
    return a_new, b_new, z_new


def pollard_rho_mod128(g, h, order, r=7, max_iterations=50000):
    
    # Random starting point: z = g^a * h^b
    a = random.randint(0, order - 1)
    b = random.randint(0, order - 1)
    
    # Init z = g^a * h^b
    g_a = (g ** a).reduce()
    h_b = (h ** b).reduce()
    z = (g_a * h_b).reduce()

    # Init points of turtoise and hare
    A_tortoise, B_tortoise, Z_tortoise = a, b, z.copy()
    A_hare, B_hare, Z_hare = a, b, z.copy()

    for i in range(max_iterations):
        # Tortoise: one step
        A_tortoise, B_tortoise, Z_tortoise = simple_update(A_tortoise, B_tortoise, Z_tortoise, g, h, order, r)
        
        # Hare: two steps
        A_hare, B_hare, Z_hare = simple_update(A_hare, B_hare, Z_hare, g, h, order, r)
        A_hare, B_hare, Z_hare = simple_update(A_hare, B_hare, Z_hare, g, h, order, r)
        
        # Check for collision
        if Z_tortoise == Z_hare:
            print(f"Collision found at iteration {i+1}")
            print(f"Tortoise: a={A_tortoise}, b={B_tortoise}")
            print(f"Hare: a={A_hare}, b={B_hare}")
            
            return
    raise Exception("No collision found within iteration limit")

# Test with smaller prime
q = 13
g = Mod128Minus3Element.from_int(2)
x_secret = 5
h = (g ** x_secret).reduce()
print(f"Testing: g = {g.to_int()}, x_secret = {x_secret}, h = {h.to_int()}")
x_found = pollard_rho_mod128(g, h, q, r=7)
print(f"x_found = {x_found}")
assert x_found == x_secret, f"Expected {x_secret}, but got {x_found}"
print("Test passed!")