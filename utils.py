def modinv(a, p):
    """Compute the modular inverse of a modulo p."""
    return pow(a, -1, p)  # Python ≥ 3.8
