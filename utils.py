def modinv(a, p):
    """Compute the modular inverse of a modulo p."""
    return pow(a, -1, p)  # Python ≥ 3.8


def modinv2(a, m):
    def egcd(a, m):
        if a == 0:
            return m, 0, 1
        g, x, y = egcd(m % a, a)
        return g, y - (m // a) * x, x
    g, x, _ = egcd(a % m, m)
    if g != 1:
        raise ValueError("Modular inverse does not exist")
    return x % m
assert modinv(3, 101) * 3 % 101 == 1  # Should pass