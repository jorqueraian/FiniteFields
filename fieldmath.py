import re

class Field:
    """
    Abstraction for fields.
    A field must satisfy a list of axioms. The axioms can be found here https://mathworld.wolfram.com/FieldAxioms.html
    In short this means for every field we must have a definition for addition, subtraction, multiplication, and division( or at least multiplicative inverses). To satisfy both subtraction and division we will define additive and multiplicative inverses. This also means we have to have a zero and an one element in our field, and these cannot be equal. The following class will act as a base class for field implementation and will be an object of operations such that for a field f. we can use the operations as follows: f.add(x, y).
    """

    def zero(self):
        """ additive element of field. x + 0 = x"""
        raise AssertionError("Not implemented")

    def one(self):
        """ multiplicative element of field. x * 1 = x"""
        raise AssertionError("Not implemented")

    def equals(self, x, y):
        """ Checks is x and y are the same element of field"""
        raise AssertionError("Not implemented")

    def add(self, x, y):
        """ adds x and y"""
        raise AssertionError("Not implemented")

    def negate(self, x):
        """ additive inverse element of field of the element x"""
        raise AssertionError("Not implemented")

    def subtract(self, x, y):
        """ subtract x and y"""
        raise AssertionError("Not implemented")

    def multiply(self, x, y):
        """ multiply x and y"""
        raise AssertionError("Not implemented")

    def reciprocal(self, x):
        """ Multiplicative inverse element of field of the element x"""
        raise AssertionError("Not implemented")

    def divide(self, x, y):
        """ multiply x and y"""
        raise AssertionError("Not implemented")
    
    def involve(self, x):
        return x

    def __eq__(self, other):
        raise AssertionError("Not implemented")

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def print_elm(self, x):
        return x


class FieldExtension(Field):
    def __init__(self, base_field, alpha_relation, involution_pow=None):
        """
            example FieldExtension(Z_5, [1,1,1])
            Z_5/(1+x+x^2)  =~ F_25

            This will be slow but will get the job done
        """
        self._bf = base_field
        self._edeg = len(alpha_relation)-1
        self._ext_rel = alpha_relation
        self.size = self._bf.size**self._edeg

        if involution_pow is not None:
            self.involve = lambda a: pow_over_field(self.is_valid(a), involution_pow, self, False)

        import polynomial as mypoly
        self._poly_helper = mypoly.PolynomialOperations(self._bf)

    def is_valid(self, x):
        if isinstance(x, str):
            x_str = x
            terms = [re.split(r'a\^?', term) for term in re.split(r'\+', x_str)]
            helper_f = lambda a: int(a[1]) if len(a)>1 and a[1] != "" else 1
            x_arr = [self._bf.zero()]*(1+max([helper_f(term) for term in terms]))  # sloppy but what ever
            for comp in re.split(r'\+', x_str):
                coeff_i = re.split(r'a\^?', comp)
                if len(coeff_i) == 1:
                    x_arr[0] = int(coeff_i[0])
                elif coeff_i[1] == "":
                    x_arr[1] = int(coeff_i[0])
                else:
                    x_arr[int(coeff_i[1])] = int(coeff_i[0])
        else:
            x_arr = list(x)
        x_arr = self._poly_helper.poly_trim(x_arr)

        for deg_monomial in range(len(x_arr)-1,self._edeg-1,-1):
            x_arr = self._poly_helper.poly_add(
                x_arr,
                [self._bf.zero()]*(deg_monomial-self._edeg)+[self._bf.multiply(self._bf.negate(x_arr[deg_monomial]), a) for a in self._ext_rel[:-1]]
            )
            x_arr[deg_monomial] = self._bf.zero()
        
        x_arr = self._poly_helper.poly_trim(x_arr)

        return x_arr

    def zero(self):
        return [self._bf.zero()]

    def one(self):
        return [self._bf.one()]

    def equals(self, x, y):
        return self._poly_helper.poly_equal(self.is_valid(x),self.is_valid(y))

    def add(self, x, y):
        return self._poly_helper.poly_add(self.is_valid(x),self.is_valid(y))

    def negate(self, x):
        x_ = self.is_valid(x)
        neg_one = self._bf.negate(self._bf.one())
        return [self._bf.multiply(neg_one, a) for a in x_]

    def subtract(self, x, y):
        return self._poly_helper.poly_subtract(self.is_valid(x), self.is_valid(y))

    def multiply(self, x, y):
        x_ = self.is_valid(x)
        y_ = self.is_valid(y)
        new_elem = [self._bf.zero()]*(self._poly_helper.poly_degree(x_)+self._poly_helper.poly_degree(y_)+1)
        for i,a in enumerate(x_):
            for j,b in enumerate(y_):
                new_elem[i+j] = self._bf.add(new_elem[i+j], self._bf.multiply(a,b))
        
        return self.is_valid(new_elem)

    def reciprocal(self, x):
        return pow_over_field(self.is_valid(x), self.size - 2, self, False)

    def divide(self, x, y):
        raise AssertionError("Not implemented")

    def __eq__(self, other):
        return True
    
    def print_elm(self, x):
        return "+".join([f"{term}{'a' if i>0 else ''}{f'^{i}' if i >1 else ''}" for i, term in enumerate(self.is_valid(x))])


class BinaryField(Field):
    """ Galois field, Binary Field. The ideas behind the implementation of each operations can be found from
    https://sites.math.washington.edu/~morrow/336_12/papers/juan.pdf
    https://en.wikipedia.org/wiki/Finite_field_arithmetic
    In General each element can be thought of as a polynomial with the ones and zeros being the coefficients"""

    def __init__(self, mod):
        """ Creates GF(2^m/mos)
        :param mod: a binary polynomial modulus. Recommended to use 0x11d
        """
        self.mod = mod
        self.size = 1 << (mod.bit_length() - 1)

        self.mults = [[None] * self.size for _ in range(self.size)]
        self.recips = [None] * self.size
        self.pows = [[None] * self.size for _ in range(self.size)]
        self._init_tables()

    def _init_tables(self):
        for i in range(0, self.size):
            for j in range(0, self.size):
                self.mults[i][j] = self.multiply(i, j)

        for i in range(self.size):
            for j in range(self.size):
                self.pows[i][j] = pow_over_field(i, j, self)

        for i in range(1, self.size):
            self.recips[i] = self.reciprocal(i)

    def is_valid(self, x):
        """ Each element of out field will be an Int and must be less than the size. Returns the element"""
        if not isinstance(x, int):
            raise TypeError()
        if x < 0 or x >= self.size:
            raise ValueError(f"{x} is not an element of the field.")
        return x

    def zero(self):
        """ additive identity element of field"""
        return 0

    def one(self):
        """ multiplicative identity element of field"""
        return 1

    def equals(self, x, y):
        """ Checks is x and y are the same element of field"""
        return self.is_valid(x) == self.is_valid(y)

    def add(self, x, y):
        """ adds x and y. Similar to adding polynomial we add coefficient. Where the coefficients live in the GF(2)
         So addition is simply for each coeff c_k = x_k + y_k % 2. so the result is x ^ y"""
        return self.is_valid(x) ^ self.is_valid(y)

    def negate(self, x):
        """ additive inverse element of field of the element x. With the definition of addition, x + x = 0. so the
        additive inverse is just its self"""
        return self.is_valid(x)

    def subtract(self, x, y):
        """ subtract x and y. With the definition of addition, x + x = 0. so subtraction is just x ^ y"""
        return self.is_valid(x) ^ self.is_valid(y)

    def multiply(self, x, y):
        """ multiply x and y. Ther is the same version as found on the wikipidea page above"""
        self.is_valid(x)
        self.is_valid(y)
        p = 0
        while y != 0:
            if y & 1 != 0:
                p ^= x
            x <<= 1
            if x >= self.size:
                x ^= self.mod
            y >>= 1
        return p

    def reciprocal(self, x):
        """ Multiplicative inverse element of field of the element x. To find an inverse we will use the extended The
        Itoh-Tsujii Algorithm. Process found here: http://cse.iitkgp.ac.in/~debdeep/osscrypto/psec/psec/pubs/thesis.pdf
        """
        return pow_over_field(x, self.size - 2, self)

    def divide(self, x, y):
        """ multiply x and y"""
        raise AssertionError("Not implemented")

    def __eq__(self, other):
        if isinstance(other, BinaryField):
            if self.mod == other.mod:
                return True
        return False


class Zp(Field):
    def __init__(self, p):
        """ Creates Z_p or Z/pZ where p is prime"""
        self.p = p
        self.size = p

    def is_valid(self, x):
        """ x must be of type np.int32 and between 0 and p. I removed the validity check here to make things faster."""
        return x

    def zero(self):
        """ additive identity element of field"""
        return 0

    def one(self):
        """ multiplicative identity element of field"""
        return 1

    def equals(self, x, y):
        """ Checks is x and y are the same element of field"""
        return self.is_valid(x) == self.is_valid(y)

    def add(self, x, y):
        return (self.is_valid(x) + self.is_valid(y)) % self.p

    def negate(self, x):
        """ additive inverse element of field of the element x """
        return (-1 * self.is_valid(x)) % self.p

    def subtract(self, x, y):
        return (self.is_valid(x) - self.is_valid(y)) % self.p

    def multiply(self, x, y):
        return (self.is_valid(x) * self.is_valid(y)) % self.p

    def reciprocal(self, x):
        """ Multiplicative inverse element of field of the element x """
        def gcd_extended(a, b):
            if a == 0:
                return b, 0, 1

            gcd, x1, y1 = gcd_extended(b % a, a)

            xp = y1 - (b // a) * x1
            yp = x1

            return gcd, xp, yp

        one, x_inv, p = gcd_extended(self.is_valid(x), self.p)
        assert one == 1, f"Extended Euclidean algorithm failed to find inv of {x}"
        return (x_inv % self.p)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.p == other.p:
                return True
        return False


def pow_over_field(base, exp, field, use_saved_vals=True):
    if use_saved_vals:
        if exp < 0:
            raise ValueError("Exp can not be negative.")
        if not isinstance(field, Field):
            raise Exception("Must provide a valid Field.")
        if isinstance(field, BinaryField) and field.pows[base][exp] is not None:
            return field.pows[base][exp]
        if isinstance(field, BinaryField) and field.pows[base][exp - 1] is not None:
            return field.multiply(base, field.pows[base][exp - 1])

    res = field.one()
    for _ in range(exp):
        res = field.multiply(base, res)
    return res


class Matrix:
    """ Matrix class that works with provided field"""
    def __init__(self, rows, columns, field, init_val=None):
        """initializes rows x columns matrix with init_val """
        if rows <= 0 or columns <= 0:
            raise ValueError("rows and columns must be > 0.")
        if not isinstance(field, Field) or not isinstance(rows, int) or not isinstance(columns, int):
            raise TypeError()

        self.f = field

        self.values = [[init_val] * columns for _ in range(rows)]
        self.rows = rows
        self.columns = columns

    def get(self, r, c):
        if r < 0 or c < 0 or r >= self.rows or c >= self.columns:
            raise IndexError("Index for r or c out of bounds")
        return self.values[r][c]

    def set(self, r, c, val):
        if r < 0 or c < 0 or r >= self.rows or c >= self.columns:
            raise IndexError("Index for r or c out of bounds")
        self.values[r][c] = self.f.is_valid(val)

    def get_sub_matrix(self, row_i, row_t, col_i, col_t):
        if row_i is None:
            row_i = 0
        if row_t is None:
            row_t = self.rows
        if col_i is None:
            col_i = 0
        if col_t is None:
            col_t = self.columns

        if row_t <= row_i or col_t <= col_i:
            raise Exception("Invalid parameters. Terminator can not be leq init.")
        result = self.__class__(row_t - row_i, col_t - col_i, self.f)
        for r in range(row_i, row_t):
            for c in range(col_i, col_t):
                result.set(r - row_i, c - col_i, self.get(r, c))
        return result

    def to_list(self, single=False):
        lst = []
        for r in range(self.rows):
            row = []
            for c in range(self.columns):
                if single:
                    lst.append(self.get(r, c))
                else:
                    row.append(self.get(r, c))
            if not single:
                lst.append(row)
        return lst

    def __str__(self):
      matrix_str = "   "
      pfn = self.f.print_elm
      spacing = 2+max([max([len(str(pfn(elm))) for elm in row]) for row in self.values])
      for (i, row) in enumerate(self.values):
          if i > 0:
              matrix_str += " \n    "
          matrix_str += " ".join(str(pfn(val)) + " "*(spacing-len(str(pfn(val)))) for val in row)
      return matrix_str + ""

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError()
        if self.f != other.f:
            raise Exception("Fields do align.")

        if self.columns != other.rows:
            raise Exception("Can not multiple matrices, inner dimensions do no align.")

        result = self.__class__(self.rows, other.columns, self.f)
        for r in range(result.rows):
            for c in range(result.columns):
                val = self.f.zero()
                for i in range(self.columns):
                    val = self.f.add(val, self.f.multiply(self.get(r, i), other.get(i, c)))
                result.set(r, c, val)
        return result


    def __add__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError()
        if self.f != other.f:
            raise Exception("Fields do align.")

        if self.columns != other.columns or self.rows != other.rows:
            raise Exception("Can not add matrices, dimensions do no align.")

        result = self.__class__(self.rows, self.columns, self.f)
        for r in range(result.rows):
            for c in range(result.columns):
                val = self.f.add(self.get(r, c), other.get(r, c))
                result.set(r, c, val)
        return result

    def __sub__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError()
        if self.f != other.f:
            raise Exception("Fields do align.")

        if self.columns != other.columns or self.rows != other.rows:
            raise Exception("Can not add matrices, dimensions do no align.")

        result = self.__class__(self.rows, self.columns, self.f)
        for r in range(result.rows):
            for c in range(result.columns):
                val = self.f.add(self.get(r, c), other.f.negate(other.get(r, c)))
                result.set(r, c, val)
        return result

    def transpose(self):
        result = self.__class__(self.columns, self.rows, self.f)
        for r in range(result.rows):
            for c in range(result.columns):
                result.set(r, c, self.get(c, r))
        return result
    
    def conj_transpose(self):
        result = self.__class__(self.columns, self.rows, self.f)
        for r in range(result.rows):
            for c in range(result.columns):
                result.set(r, c, self.f.involve(self.get(c, r)))
        return result

    def any(self):
        for r in range(self.rows):
            for c in range(self.columns):
                if not self.f.equals(self.get(r, c), self.f.zero()):
                    return True
        return False

    def copy(self):
        result = self.__class__(self.rows, self.columns, self.f)
        result.values = self.values
        return result

    def swap_rows(self, r1, r2):
        if r1 < 0 or r2 < 0 or r1 >= self.rows or r2 >= self.rows:
            raise IndexError("Index for r or c out of bounds")
        self.values[r1], self.values[r2] = self.values[r2], self.values[r1]

    def multiply_row(self, r, multiplier):
        if r < 0 or r >= self.rows:
            raise IndexError("Index for r or c out of bounds")
        self.values[r] = [self.f.multiply(multiplier, val) for val in self.values[r]]

    def add_rows(self, r1, r2, multiplier):
        """ r2 = r2 + r1 * multiplier """
        if r1 < 0 or r2 < 0 or r1 >= self.rows or r2 >= self.rows:
            raise IndexError("Index for r or c out of bounds")
        self.values[r2] = [self.f.add(val_r2, self.f.multiply(multiplier, val_r1)) for val_r1, val_r2 in
                           zip(self.values[r1], self.values[r2])]

    def row_echelon_form(self):
        lead = 0
        for r in range(self.rows):
            if lead >= self.columns:
                return
            i = r
            while self.values[i][lead] == 0:
                i += 1
                if i == self.rows:
                    i = r
                    lead += 1
                    if self.columns == lead:
                        return
            self.swap_rows(i, r)
            lv_recip = self.f.reciprocal(self.values[r][lead])

            for i in range(r, self.rows):
                if i != r:
                    lv = self.values[i][lead]
                    self.add_rows(r, i, self.f.negate(self.f.multiply(lv, lv_recip)))
            lead += 1

    def reduced_row_echelon_form(self):
        lead = 0
        for r in range(self.rows):
            if lead >= self.columns:
                return
            i = r
            while self.values[i][lead] == self.f.zero():
                i += 1
                if i == self.rows:
                    i = r
                    lead += 1
                    if self.columns == lead:
                        return
            self.swap_rows(i, r)
            lv_recip = self.f.reciprocal(self.values[r][lead])
            self.multiply_row(r, lv_recip)

            for i in range(self.rows):
                if i != r:
                    lv = self.values[i][lead]
                    self.add_rows(r, i, self.f.negate(lv))
            lead += 1

    def kernel_space(self):
        """
        I used a rather simple method for doing this as explained in the Computation by Gaussian elimination section
        found here https://en.wikipedia.org/wiki/Kernel_(linear_algebra).
        """
        ai_matrix = augmented_a_b_matrix(self.transpose(), identity_n(self.columns, self.f))
        ai_matrix.reduced_row_echelon_form()
        res = []
        for r in range(ai_matrix.rows):
            valid = True
            for c in range(self.rows):
                if not self.f.equals(ai_matrix.get(r, c), self.f.zero()):
                    valid = False
                    break
            if valid:
                res.append([ai_matrix.get(r, c) for c in range(self.rows, ai_matrix.columns)])

        if len(res) == 0 or len(res[0]) == 0:
            return 0
        result = self.__class__(len(res), len(res[0]), self.f)
        for r in range(len(res)):
            for c in range(len(res[0])):
                result.set(r, c, res[r][c])
        return result.transpose()


def identity_n(n, field, val=None):
    """
    returns a nxn identity matrix
    """
    if not isinstance(field, Field):
        raise TypeError()

    result = Matrix(n, n, field)
    for r in range(result.rows):
        for c in range(result.columns):
            if r == c:
                if val is None:
                    result.set(r, c, field.one())
                else:
                    result.set(r, c, val)
            else:
                result.set(r, c, field.zero())
    return result


def augmented_a_b_matrix(a, b):
    """
    Simply combines to matrices as a augmented matrix. Returns (A | B)
    """
    if not isinstance(a, Matrix) or not isinstance(b, Matrix):
        raise TypeError()
    if a.f != b.f:
        raise Exception("Fields do no align.")

    if a.rows != b.rows:
        raise Exception(f"Matrices do no align. a: {a.rows}x{a.columns} and b: {b.rows}x{b.columns}")

    axb = Matrix(a.rows, a.columns + b.columns, a.f)
    for r in range(axb.rows):
        for c in range(axb.columns):
            if c >= a.columns:
                val = b.get(r, c - a.columns)
            else:
                val = a.get(r, c)
            axb.set(r, c, val)
    return axb


def solve_ax_b(a, b):
    """
    Creates an augmented matrix and find the RREF. After that parses solution. Note this does not check for singular
    matrices and tries to provide a particular solution, may have unexpected behavior.
    """

    if not isinstance(a, Matrix) or not isinstance(b, Matrix):
        raise TypeError
    if b.columns != 1 or b.rows != a.rows:
        raise Exception("Matrix b must be nx1.")

    axb = augmented_a_b_matrix(a, b)
    axb.row_echelon_form()

    # now we got to find the solution
    c = 0
    result = Matrix(a.columns, 1, a.f, init_val=a.f.zero())
    if axb.rows < axb.columns - 1:
        raise Exception("Can not solve, too few rows")
    for r in range(axb.rows-1, -1, -1):
        c = (min(axb.columns - 2, r))
        if c != r:
            if axb.get(r, c) != 0 or axb.get(r, c + 1) != 0:
                raise Exception("Inconsistent linear system, or singular matrix A")
        else:
            to_sub = axb.f.zero()
            for c2 in range(c + 1, axb.columns - 1):
                to_sub = axb.f.add(to_sub, axb.f.multiply(axb.get(r, c2), result.get(c2, 0)))
            result.set(c, 0, axb.f.multiply(axb.f.reciprocal(axb.get(r, c)),
                                            axb.f.subtract(axb.get(r, axb.columns - 1), to_sub)))
    return result


def solve_lstsq(a, b):
    """
    This function solves the normal equations A^T*A*x = A^T*b. Not this function considers only the 'euclidean dot
    product'
    This really doesnt work with the Binary field. In fact this doesnt really work at all for fields
    """
    if not isinstance(a, Matrix) or not isinstance(a, Matrix):
        raise TypeError

    ata = a.transpose() * a
    atb = a.transpose() * b

    x = solve_ax_b(ata, atb)
    return x


def create_matrix(lst, field):
    """
    Helper function to more easily initialize a matrix from a list
    """
    rows = len(lst)
    if rows == 0:
        raise Exception("Invalid Input")
    if isinstance(lst[0], list):
        columns = len(lst[0])
    else:
        columns = 1

    result = Matrix(rows, columns, field)
    for r in range(rows):
        for c in range(columns):
            result.set(r, c, lst[r][c])
    return result



