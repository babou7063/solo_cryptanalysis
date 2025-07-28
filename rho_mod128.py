from mod128_minus3 import Mod128Minus3Element
import random
from utils import modinv2

def partition(z, r=7):
    return z.to_int() % r

def update(a, b, z, g, h, order, r=7):
    zone = partition(z, r)
    if zone == 0:
        a_new = (a + 1) % order
        b_new = b
        z_new = (z * g).reduce()
    elif zone == 1:
        a_new = a
        b_new = (b + 1) % order
        z_new = (z * h).reduce()
    elif zone == 2:
        a_new = (2 * a) % order
        b_new = (2 * b) % order
        z_new = z.square().reduce()
    elif zone == 3:
        a_new = (a + 2) % order
        b_new = b
        z_new = (z * g * g).reduce()
    elif zone == 4:
        a_new = a
        b_new = (b + 2) % order
        z_new = (z * h * h).reduce()
    elif zone == 5:
        a_new = (3 * a) % order
        b_new = (3 * b) % order
        z_new = (z * z * z).reduce()
    else:
        a_new = (a + b) % order
        b_new = (b + a) % order
        z_new = (z * g * h).reduce()
    return a_new, b_new, z_new

def pollard_rho_mod128(g, h, order, r=7):
    a = random.randint(0, order - 1)
    b = random.randint(0, order - 1)
    z = (g ** a * h ** b).reduce()

    A, B, Z = a, b, z.copy()
    A2, B2, Z2 = a, b, z.copy()

    for i in range(1000):  # Smaller iteration limit for small q
        A, B, Z = update(A, B, Z, g, h, order, r)
        for _ in range(2):
            A2, B2, Z2 = update(A2, B2, Z2, g, h, order, r)
        if Z.to_int() == Z2.to_int():
            if B == B2:
                print("Failure: b - b' ≡ 0, retrying")
                return pollard_rho_mod128(g, h, order, r)
            x = ((A - A2) * modinv2(B2 - B, order)) % order
            return x
    raise Exception("No collision found within iteration limit")

# Test with smaller prime
q = 101  # Small prime for quick testing
g = Mod128Minus3Element.from_int(7)
x_secret = random.randint(1, q - 1)
h = (g ** x_secret).reduce()
print(f"Testing with q = {q}, g = 7, x_secret = {x_secret}")
x_found = pollard_rho_mod128(g, h, q, r=7)
print(f"x_found = {x_found}")
assert x_found == x_secret, f"Expected {x_secret}, but got {x_found}"
print("Test passed!")