import fieldmath


# Polynomial operations. ALl functions work under some field and all parameters are
# polynomials such that 3 + 3x + x^4 is represented with the array [3, 3, 0, 0, 1]
class PolynomialOperations(object):
    def __init__(self, field):
        if isinstance(field, fieldmath.Field):
            self.f = field
        else:
            raise Exception("Not a valid field")

    def poly_call(self, poly, x):
        res = self.f.zero()
        for c in range(len(poly)-1, -1, -1):
            res = self.f.add(self.f.multiply(res, x), poly[c])
        return res

    def poly_degree(self, poly):
        """
        Return the degree of the polynomial. ex 3 + 3x + x^4 or [3, 3, 0, 0, 1] has degree 4
        """
        poly = self.poly_trim(poly)
        if poly:
            return len(poly) - 1
        else:
            return 0

    def poly_trim(self, poly):
        """
        removes all leading 0 terms from the coefficients array
        """
        while poly and self.f.equals(poly[-1], self.f.zero()):
            poly.pop()
        if len(poly) == 0:
            poly = [self.f.zero()]
        return poly

    def poly_lead(self, poly):
        """
        Returns the leading term. ex 3 + 3x + x^4 or [3, 3, 0, 0, 1] will return 1.
        """
        poly = self.poly_trim(poly)
        for i in range(len(poly) - 1, -1, -1):
            if poly[i] != 0:
                return poly[i]
        return self.f.zero()
    
    def set(self, poly, degree_term, coeff_val):
        """
        Sets coefficient of the degree term to the coeff value. ex set(2, 2) for x+2 will result in the polynomial being
        2x^2 + x + 2
        """
        if len(poly) <= degree_term:
            poly += [self.f.zero()] * (degree_term - len(poly) + 1)
        poly[degree_term] = coeff_val

    def poly_divmod(self, poly1, poly2):
        """
        Heavily modified version of a polynomial long division originally found here:
        https://rosettacode.org/wiki/Polynomial_long_division
        """
        # top = self.__class__(*self.coefficients, f=self.f, flipped=True)
        # bottom = self.__class__(*other.coefficients, f=other.f, flipped=True)
        deg_top = self.poly_degree(poly1)
        deg_bot = self.poly_degree(poly2)

        if deg_bot < 0:
            raise ZeroDivisionError

        if deg_top >= deg_bot:
            q = [self.f.zero()] * deg_top
            while deg_top >= deg_bot and self.poly_lead(poly1) != self.f.zero():  # and top.lead() != self.f.zero():
                d = [self.f.zero()] * (deg_top - deg_bot) + poly2
                inv = self.f.reciprocal(self.poly_lead(d))
                multiplier = self.f.multiply(self.poly_lead(poly1), inv)
                self.set(q, deg_top - deg_bot, multiplier)
                d = [self.f.multiply(multiplier, c) for c in d]
                poly1 = self.poly_subtract(poly1, d)
                deg_top = self.poly_degree(poly1)
        else:
            q = [self.f.zero()]
        return q, poly1

    def poly_add(self, poly1, poly2):
        return [self.f.add(t1, t2) for t1, t2 in zip_longest(poly1, poly2, fillvalue=self.f.zero())]

    def poly_subtract(self, poly1, poly2):
        return [self.f.subtract(t1, t2) for t1, t2 in zip_longest(poly1, poly2, fillvalue=self.f.zero())]

    def poly_equal(self, poly1, poly2):
        poly1 = self.poly_trim(poly1)
        poly2 = self.poly_trim(poly2)
        if len(poly1) == len(poly2):
            for e1, e2 in zip(poly1, poly2):
                if not self.f.equals(e1, e2):
                    return False
            return True
        else:
            return False

    def poly_not_equal(self, poly1, poly2):
        return not self.poly_equal(poly1, poly2)

    def poly_gcd(self, poly1, poly2):
        if self.poly_degree(poly1) > self.poly_degree(poly2):
            return self.poly_gcd(poly2, poly1)

        if self.poly_equal(poly1, [self.f.zero()]):
            return poly2

        q, r = self.poly_divmod(poly2, poly1)
        r = self.poly_trim(r)

        return self.poly_gcd(r, poly1)
    
    def poly_derivative(self, poly):
        ret_pol = [0]*self.poly_degree(poly)
        for i in range(len(ret_pol)):
            ret_pol[i] = self.f.multiply(poly[i+1], (self.f.add(i, self.f.one())))
        return ret_pol


def zip_longest(iter1, iter2, fillvalue=None):
    # I did not create this. Found with original framework for Polynomial class. used for __add__ and __sub__
    for i in range(max(len(iter1), len(iter2))):
        if i >= len(iter1):
            yield (fillvalue, iter2[i])
        elif i >= len(iter2):
            yield (iter1[i], fillvalue)
        else:
            yield (iter1[i], iter2[i])
        i += 1