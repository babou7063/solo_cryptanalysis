

# Définir une courbe elliptique simplifiée
# y² = x³ + ax + b sur ℤ/pℤ
# p = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31




# Définir les points P et Q
# P est le générateur
# Q = k·P pour un k secret (que tu caches)



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
        
        # Case of an infinity point
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
    
    def scalar_mul(self, p, k):
        pass



# Test 
E = EllipticCurve(a=2, b=3, p=97)
P = Point(3, 6, E)
O = Point.infinity(E)

print(P)
print(-P)
print(P + O)  # doit afficher (3, 6)
print(P + (-P))  # doit afficher ∞
print(P + P)