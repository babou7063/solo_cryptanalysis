"""
Defines classes for working with elliptic curves over finite fields (modulo a prime p).

Includes:
- class Point: represents a point (x, y) on the curve or the point at infinity
- class EllipticCurve: represents a short Weierstrass elliptic curve of the form y² = x³ + ax + b over ℤ/pℤ, and provides point addition and scalar multiplication
"""


class Point:
    def __init__(self, x, y, curve, at_infinity=False):
        """Create a new point on a curve.
        
        :param x: x-coordinate
        :param y: y-coordinate
        :param curve: the elliptic curve to create the point on
        :param at_infinity: whether the point is the infinity point
        """
        
        self.x = x
        self.y = y
        self.curve = curve
        self.at_infinity = at_infinity

    @staticmethod
    def infinity(curve):
        """Return the point at infinity on the given curve.
        
        :param curve: the elliptic curve to create the point on
        :return: the point at infinity
        """
        return Point(None, None, curve, at_infinity=True)

    def __eq__(self, other):
        """Compare two points for equality."""
        
        if self.at_infinity and other.at_infinity:
            return True
        if self.at_infinity or other.at_infinity:
            return False
        return self.x == other.x and self.y == other.y

    def __neg__(self):
        """Return the negation of the point.
        
        This is the point on the same x-coordinate but with the negative of the y-coordinate.
        """
        if self.at_infinity:
            return self
        return Point(self.x, (-self.y) % self.curve.p, self.curve)

    def __add__(self, other):
        """Add two points on the same curve. """

        return self.curve.add_points(self, other)

    def __str__(self):
        """Return a string representation of the point."""
        
        if self.at_infinity:
            return "Inf"
        return f"({self.x}, {self.y})"
    



class EllipticCurve:
    
    def __init__(self, a, b, p):
        """Initialize an elliptic curve given by the equation y^2 = x^3 + ax + b over the finite field of integers modulo p.

        :param a: coefficient a in the curve equation
        :param b: coefficient b in the curve equation
        :param p: prime number specifying the finite field
        """

        self.a, self.b, self.p = a, b, p
    
    def add_points(self, P, Q):
        """Addition of two points on the elliptic curve.

        :param P: first point
        :param Q: second point
        :return: resulting point
        """
        
        # Case of an infinity point
        if P.at_infinity and Q.at_infinity:
            return Point.infinity(self)
        if P.at_infinity:
            return Q
        if Q.at_infinity:
            return P

        # Case of P + (-P) or P = Q but y = 0
        if P.x == Q.x and (P.y != Q.y or P.y == 0):
            return Point.infinity(self)

        if P == Q:
            num = (3 * P.x * P.x + self.a) % self.p   # numerator
            den = pow(2 * P.y, -1, self.p)            # modulus inverse of denominator
        else:
            # Classic addition
            num = (Q.y - P.y) % self.p                # numerator
            den = pow(Q.x - P.x, -1, self.p)          # modulus inverse of denominator

        # Compute the coordinates of the resulting point
        lam = (num * den) % self.p
        x_r = (lam * lam - P.x - Q.x) % self.p
        y_r = (lam * (P.x - x_r) - P.y) % self.p

        return Point(x_r, y_r, self)
    
    def scalar_mul(self, k, P):
        """Perform scalar multiplication of a point on 
        the elliptic curve with the double-and-add method

        :param k: integer scalar to multiply
        :param P: Point on the elliptic curve
        :return: Resulting point after multiplying P by k
        """

        if k % self.p == 0 or P.at_infinity:
            return Point.infinity(self)

        result = Point.infinity(self)
        addend = P

        while k > 0:
            if k & 1:  # if bit is 1, add the addend to the result if k is odd
                result = result + addend
            addend = addend + addend
            k >>= 1  # divide by 2

        return result
    
    def find_order(self, P, max_iter=1000):
        """Compute the order of a point on the elliptic curve.
        
        :param P: point on the elliptic curve
        :param max_iter: maximum number of iterations to compute the order (default: 1000), limite to avoid loops
        :return: the order of the point
        :raises ValueError: if the order is not found within max_iter iterations
        """
        assert not P.at_infinity, "Cannot compute order of the point at infinity."
        
        R = P
        for n in range(1, max_iter + 1):
            if R.at_infinity:
                return n
            R = R + P
        raise ValueError("Point order not found within max_iter")




""" Test
p = 97
E = EllipticCurve(a=2, b=3, p=p)
P = Point(3, 9, E)

k = 5
Q = E.scalar_mul(k, P)
print(f"{k} * P = {Q}")


Q_manual = P + P + P + P + P
print(f"Manual check : {Q_manual}")

p = 97
a = 2
b = 3
curve = EllipticCurve(a, b, p)
P = Point(95, 6, curve)
order = curve.find_order(P)
print(f"Order of P: {order}")
"""
