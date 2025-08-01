
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
        """Convert an integer to a Mod128Minus3Element.
        Uses base-2^12.8 decomposition.
        
        :param n: the integer to convert
        :return: the corresponding Mod128Minus3Element
        """
        coeffs = []
        for exp in Mod128Minus3Element.BASE_EXPONENTS:
            mask = (1 << 13) - 1
            coeffs.append(n & mask)
            n >>= 13
        return Mod128Minus3Element(coeffs)

    def to_int(self): 
        """Reconstruct the integer value (not reduced modulo 2^128 - 3).
        
        :return: Integer representation of the element
        """

        result = 0
        for i, fi in enumerate(self.coeffs):
            result += fi << self.BASE_EXPONENTS[i]
        return result
    
    def __eq__(self, other):
        return self.reduce().to_int() == other.reduce().to_int()

    def __add__(self, other):
        """Add two elements coefficient-wise.

        :param other: the other element to add
        :return: the sum of the two elements
        """
        return Mod128Minus3Element([
            self.coeffs[i] + other.coeffs[i] for i in range(10)
        ])

    def __sub__(self, other):
        """Subtract two elements coefficient-wise.

        :param other: the other element to subtract
        :return: the difference of the two elements
        """
        return Mod128Minus3Element([
            self.coeffs[i] - other.coeffs[i] for i in range(10)
        ])

    # Multiplication algorithm from Bernstein, Lange, Schwabe, 2010
    # See: https://dl.acm.org/doi/10.5555/1964658.1964669
    def __mul__(self, other):    
        """
        Multiply two elements modulo 2^128 - 3 (before reduction).
        Returns a new Mod128Minus3Element instance.

        :param other: The Mod128Minus3Element to multiply with
        :return: A new Mod128Minus3Element representing the product
        """

        f = self.coeffs
        g = other.coeffs
        r = [0] * 10  # result

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
    
    # Reduction strategy inspired by the carry-chain optimization in the paper of Bernstein, Lange, Schwabe, 2010
    # See: https://dl.acm.org/doi/10.5555/1964658.1964669
    def reduce(self):
        """
        TODO : name from where this method comes
        Reduce coefficients using carry chain from r0 -> r1 -> ... -> r9 -> r0 -> r1
        Ensures all coefficients are bounded appropriately.

        :return: this element, with coefficients reduced modulo 2^128 - 3
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

        :return: A new Mod128Minus3Element representing the square of self
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
    
    def __pow__(self, exponent, modulo=None):
        """
        Modular exponentiation (right-to-left binary method).

        :param exponent: The exponent to raise the element to
        :param modulo: Ignored, here for compatibility with built-in pow()
        :return: The result of self raised to the power of exponent, modulo 2^128 - 3
        """
        if exponent < 0:
            raise ValueError("Negative exponents not supported (modular inverse not defined here).")

        result = Mod128Minus3Element.from_int(1)
        base = self.copy()

        while exponent > 0:
            if exponent % 2 == 1:
                result = (result * base).reduce()
            base = base.square()
            exponent >>= 1

        return result

class Point:
    def __init__(self, x: Mod128Minus3Element, y: Mod128Minus3Element, infinity=False):
        self.x = x
        self.y = y
        self.infinity = infinity
    
    @classmethod
    def at_infinity(cls):
        return Point(Mod128Minus3Element.zero(), Mod128Minus3Element.zero(), infinity=True)
    
    def __eq__(self, other):
        if self.at_infinity and other.at_infinity:
            return True
        if self.at_infinity or other.at_infinity:
            return False
        return self.x == other.x and self.y == other.y
    
    def copy(self):
        return Point(self.x.copy(), self.y.copy())
    
    def negate(self):
        """Return -P = (x, -y)"""
        return Point(self.x.copy(), Mod128Minus3Element.zero() - self.y)
    
    def canonical_form(self):
        """Return canonical form |P| where we choose between P and -P"""
        # Use y-coordinate parity as in the paper
        y_int = self.y.to_int()
        if y_int % 2 == 0:
            return self.copy()
        else:
            return self.negate()
    
    def double(self):
        """
        Double the point using the affine coordinates formula for y² = x³ + ax + b:
        For the curve y² = x³ - x (a = -1, b = 0):
        λ = (3x² - 1)/(2y)
        x3 = λ² - 2x
        y3 = λ(x - x3) - y
        
        :return: a new Point, being the double of self
        :raises ValueError: if the denominator of the double formula is zero
        """
        # Case of point at infinity
        if self.infinity:
            return Point.at_infinity()
        
        # Check if point is (0, 0) - this would be at infinity for most curves
        zero = Mod128Minus3Element.from_int(0)
        if self.x == zero and self.y == zero:
            return Point.at_infinity()
        
        # Check if y = 0 (point is its own inverse, so 2P = O)
        if self.y == zero:
            return Point.at_infinity()
        
        # Compute slope λ = (3x² - 1)/(2y)
        # Numerator: 3x² - 1
        x_squared = self.x.square()
        three_x_squared = x_squared * Mod128Minus3Element.from_int(3)
        numerator = three_x_squared - Mod128Minus3Element.from_int(1)
        
        # Denominator: 2y
        denominator = self.y * Mod128Minus3Element.from_int(2)
        
        # Compute modular inverse of denominator
        p = (2**128 - 3) // 76439
        denom_int = denominator.to_int() % p
        
        if denom_int == 0:
            return Point.at_infinity()
        
        # Compute λ = numerator * denominator^(-1)
        denom_inv = pow(denom_int, p - 2, p)
        lambda_int = (numerator.to_int() * denom_inv) % p
        lambda_elem = Mod128Minus3Element.from_int(lambda_int)
        
        # x3 = λ² - 2x
        lambda_squared = lambda_elem.square()
        two_x = self.x * Mod128Minus3Element.from_int(2)
        x3 = lambda_squared - two_x
        
        # y3 = λ(x - x3) - y
        x_diff = self.x - x3
        y3 = (lambda_elem * x_diff) - self.y
        
        return Point(x3.reduce(), y3.reduce())


    def add(self, other):
        """Add two points on the elliptic curve"""
        # Handle infinity cases
        if self.infinity and other.infinity:
            return Point.at_infinity()
        if self.infinity:
            return other.copy()
        if other.infinity:
            return self.copy()
        
        # Check if points are equal
        if self == other:
            return self.double()
        
        # Check if points are inverses (same x, opposite y)
        if self.x == other.x:
            if self.y == (Mod128Minus3Element.zero() - other.y):
                return Point.at_infinity()
            # If same x but different y (and not opposites), this shouldn't happen on a valid curve
            return Point.at_infinity()
        
        # Compute slope λ = (y2 - y1)/(x2 - x1)
        numerator = other.y - self.y
        denominator = other.x - self.x
        
        # Compute modular inverse
        p = (2**128 - 3) // 76439
        denom_int = denominator.to_int() % p
        
        if denom_int == 0:
            return Point.at_infinity()
        
        denom_inv = pow(denom_int, p - 2, p)
        lambda_int = (numerator.to_int() * denom_inv) % p
        lambda_elem = Mod128Minus3Element.from_int(lambda_int)
        
        # x3 = λ² - x1 - x2
        lambda_squared = lambda_elem.square()
        x3 = lambda_squared - self.x - other.x
        
        # y3 = λ(x1 - x3) - y1
        x_diff = self.x - x3
        y3 = (lambda_elem * x_diff) - self.y
        
        return Point(x3.reduce(), y3.reduce())




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

# Double
P = Point(Mod128Minus3Element.from_int(12345), Mod128Minus3Element.from_int(67890))
P2 = P.double()
print(f"P = ({P.x.to_int()}, {P.y.to_int()})")
print(f"P2 = ({P2.x.to_int()}, {P2.y.to_int()})")

# Addition
P = Point(Mod128Minus3Element.from_int(12345), Mod128Minus3Element.from_int(67890))
P_add = P.add(P)
print(f"P = ({P.x.to_int()}, {P.y.to_int()})")
print(f"P + P = ({P_add.x.to_int()}, {P_add.y.to_int()})")
"""