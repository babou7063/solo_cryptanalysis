
class Mod128Minus3Element:
    """
    Represents an integer modulo 2^128 - 3 using 10-coefficient representation.
    Supports addition, multiplication, and reduction.
    """
    BASE_EXPONENTS = [0, 13, 26, 39, 52, 64, 77, 90, 103, 116]  # exponents for radix 2^12.8
    BASE = 2 ** 128 - 3
    WORD_BITS = 13  # approx base-2^12.8
    
    def __init__(self, coeffs):
        assert len(coeffs) == 10
        self.coeffs = coeffs[:]  # Copy the list

    def __repr__(self):
        return f"Mod128Minus3({self.coeffs})"
    
    def copy(self):
        return Mod128Minus3Element(self.coeffs[:])

    @staticmethod
    def zero():
        return Mod128Minus3Element([0] * 10)

    @staticmethod
    def from_int(n):
        """
        Create a Mod128Minus3Element from an integer n.
        Uses base-2^12.8 decomposition.
        """
        
        coeffs = []
        for exp in Mod128Minus3Element.BASE_EXPONENTS:
            mask = (1 << 13) - 1
            coeffs.append(n & mask)
            n >>= 13
        return Mod128Minus3Element(coeffs)

    def to_int(self):
        """
        Reconstruct the integer value (not reduced modulo 2^128 - 3).
        """
        result = 0
        for i, fi in enumerate(self.coeffs):
            result += fi << self.BASE_EXPONENTS[i]
        return result

    def __add__(self, other):
        """
        Add two elements coefficient-wise.
        """
        return Mod128Minus3Element([
            self.coeffs[i] + other.coeffs[i] for i in range(10)
        ])

    def __sub__(self, other):
        """
        Subtract two elements coefficient-wise.
        """
        return Mod128Minus3Element([
            self.coeffs[i] - other.coeffs[i] for i in range(10)
        ])
    
    def __mul__(self, other):
        """
        Multiply two elements modulo 2^128 - 3 (before reduction).
        Returns a new Mod128Minus3Element instance.
        """
        f = self.coeffs
        g = other.coeffs
        r = [0] * 10  # result

        # Multiplication algorithm from .... (mettre lien du papier ici)
        r[0] = f[0]*g[0] + 3*(2*f[9]*g[1] + 2*f[8]*g[2] + 2*f[7]*g[3] + 2*f[6]*g[4] + f[5]*g[5] + 2*f[4]*g[6] + 2*f[3]*g[7] + 2*f[2]*g[8] + 2*f[1]*g[9])
        r[1] = f[1]*g[0] + f[0]*g[1] + 3*(2*f[9]*g[2] + 2*f[8]*g[3] + 2*f[7]*g[4] + f[6]*g[5] + f[5]*g[6] + 2*f[4]*g[7] + 2*f[3]*g[8] + 2*f[2]*g[9])
        r[2] = f[2]*g[0] + f[1]*g[1] + f[0]*g[2] + 3*(2*f[9]*g[3] + 2*f[8]*g[4] + f[7]*g[5] + f[6]*g[6] + f[5]*g[7] + 2*f[4]*g[8] + 2*f[3]*g[9])
        r[3] = f[3]*g[0] + f[2]*g[1] + f[1]*g[2] + f[0]*g[3] + 3*(2*f[9]*g[4] + f[8]*g[5] + f[7]*g[6] + f[6]*g[7] + f[5]*g[8] + 2*f[4]*g[9])
        r[4] = f[4]*g[0] + f[3]*g[1] + f[2]*g[2] + f[1]*g[3] + f[0]*g[4] + 3*(f[9]*g[5] + f[8]*g[6] + f[7]*g[7] + f[6]*g[8] + f[5]*g[9])
        r[5] = f[5]*g[0] + 2*f[4]*g[1] + 2*f[3]*g[2] + 2*f[2]*g[3] + 2*f[1]*g[4] + f[0]*g[5] + 3*(2*f[9]*g[6] + 2*f[8]*g[7] + 2*f[7]*g[8] + 2*f[6]*g[9])
        r[6] = f[6]*g[0] + f[5]*g[1] + 2*f[4]*g[2] + 2*f[3]*g[3] + 2*f[2]*g[4] + f[1]*g[5] + f[0]*g[6] + 3*(2*f[9]*g[7] + 2*f[8]*g[8] + 2*f[7]*g[9])
        r[7] = f[7]*g[0] + f[6]*g[1] + f[5]*g[2] + 2*f[4]*g[3] + 2*f[3]*g[4] + f[2]*g[5] + f[1]*g[6] + f[0]*g[7] + 3*(2*f[9]*g[8] + 2*f[8]*g[9])
        r[8] = f[8]*g[0] + f[7]*g[1] + f[6]*g[2] + f[5]*g[3] + 2*f[4]*g[4] + f[3]*g[5] + f[2]*g[6] + f[1]*g[7] + f[0]*g[8] + 3*(2*f[9]*g[9])
        r[9] = f[9]*g[0] + f[8]*g[1] + f[7]*g[2] + f[6]*g[3] + f[5]*g[4] + f[4]*g[5] + f[3]*g[6] + f[2]*g[7] + f[1]*g[8] + f[0]*g[9]

        return Mod128Minus3Element(r)
    
    def reduce(self):
        """
        Reduce coefficients using carry chain from r0 -> r1 -> ... -> r9 -> r0 -> r1
        Ensures all coefficients are bounded appropriately.
        """
        r = self.coeffs[:]  # copy to avoid in-place issues

        # Carry steps r0 -> r1, ..., r3 -> r4
        for i in [0, 1, 2, 3]:
            c = (r[i] + (1 << 12)) >> 13
            r[i] -= c << 13
            r[i+1] += c

        # Special carry for r4 -> r5
        c = (r[4] + (1 << 11)) >> 12
        r[4] -= c << 12
        r[5] += c

        # Carry steps r5 -> r6 -> ... -> r8 -> r9
        for i in [5, 6, 7, 8]:
            c = (r[i] + (1 << 12)) >> 13
            r[i] -= c << 13
            r[i+1] += c

        # Final step: r9 -> r0 (modulo 2^128 - 3)
        c = (r[9] + (1 << 11)) >> 12
        r[9] -= c << 12
        r[0] += 3 * c  # Important: multiply by 3 (mod 2^128 ≡ 3)

        # Redo r0 -> r1 -> r2
        for i in [0, 1]:
            c = (r[i] + (1 << 12)) >> 13
            r[i] -= c << 13
            r[i+1] += c

        self.coeffs = r
        return self
    
    def square(self):
        """
        Optimized squaring: avoids redundant terms and reuses symmetric products.
        """
        f = self.coeffs
        r = [0] * 10

        # Pre-compute to reduce the cost
        f2 = [2 * fi for fi in f]
        f3 = [3 * fi for fi in f]
        f6 = [6 * fi for fi in f]

        # r0 to r9 (with simplification for f_i * f_j = 2*f_i*f_j si i != j)
        r[0] = f[0] * f[0] + 3 * (2 * f[9] * f[1] + 2 * f[8] * f[2] + 2 * f[7] * f[3] + 2 * f[6] * f[4] + f[5] * f[5] + 2 * f[4] * f[6] + 2 * f[3] * f[7] + 2 * f[2] * f[8] + 2 * f[1] * f[9])
        r[1] = f2[0]*f[1] + 3 * (2 * f[9] * f[2] + 2 * f[8] * f[3] + 2 * f[7] * f[4] + f[6] * f[5] + f[5] * f[6] + 2 * f[4] * f[7] + 2 * f[3] * f[8] + 2 * f[2] * f[9])
        r[2] = f2[0]*f[2] + f[1]*f[1] + 3 * (2 * f[9] * f[3] + 2 * f[8] * f[4] + f[7] * f[5] + f[6] * f[6] + f[5] * f[7] + 2 * f[4] * f[8] + 2 * f[3] * f[9])
        r[3] = f2[0]*f[3] + f2[1]*f[2] + 3 * (2 * f[9] * f[4] + f[8] * f[5] + f[7] * f[6] + f[6] * f[7] + f[5] * f[8] + 2 * f[4] * f[9])
        r[4] = f2[0]*f[4] + f2[1]*f[3] + f[2]*f[2] + 3 * (f[9] * f[5] + f[8] * f[6] + f[7] * f[7] + f[6] * f[8] + f[5] * f[9])
        r[5] = f2[0]*f[5] + f2[1]*f[4] + f2[2]*f[3] + 3 * (2 * f[9] * f[6] + 2 * f[8] * f[7] + 2 * f[7] * f[8] + 2 * f[6] * f[9])
        r[6] = f2[0]*f[6] + f2[1]*f[5] + f2[2]*f[4] + f[3]*f[3] + 3 * (2 * f[9] * f[7] + 2 * f[8] * f[8] + 2 * f[7] * f[9])
        r[7] = f2[0]*f[7] + f2[1]*f[6] + f2[2]*f[5] + f2[3]*f[4] + 3 * (2 * f[9] * f[8] + 2 * f[8] * f[9])
        r[8] = f2[0]*f[8] + f2[1]*f[7] + f2[2]*f[6] + f2[3]*f[5] + f[4]*f[4] + 3 * (2 * f[9] * f[9])
        r[9] = f2[0]*f[9] + f2[1]*f[8] + f2[2]*f[7] + f2[3]*f[6] + f2[4]*f[5]

        result = Mod128Minus3Element(r)
        return result.reduce()



######################################
#             Test                   #
######################################
"""
# Multiplication
a = Mod128Minus3Element.from_int(42)
b = Mod128Minus3Element.from_int(17)
c = a * b

print("a =", a)
print("b =", b)
print("a * b =", c)
print("int(a * b) =", c.to_int())
print("42 * 17 mod (2^128 - 3) =", (42 * 17) % (2**128 - 3))

# Reduce
a = Mod128Minus3Element([2**12, 0, 0, 0, 0, 0, 0, 0, 0, 0])
print("Before reduction:", a)
a.reduce()
print("After reduction:", a)

# Square
a = Mod128Minus3Element.from_int(123456789)
sq = a.square()
print("a =", a)
print("a^2 =", sq)
print("a^2 (int) =", sq.to_int())
print("(a.to_int() ** 2) % (2**128 - 3) =", pow(a.to_int(), 2, 2**128 - 3))
"""