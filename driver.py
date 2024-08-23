# For field implementation look at fieldmath.py
import fieldmath
import polynomial
"""# Finite field stuff. Originally created for error correcting code implementations. Revisions needed 

# This file will highlight what fieldmath.py can do
# First we will look at the Binary field of {0, 1} used for the more traditional error-correcting codes
bf = fieldmath.Zp(5)
print("Working in the finite field Z/pZ.")

t = 3

x = bf.multiply(bf.add(bf.multiply(t,t), 1), bf.reciprocal(bf.subtract(1, bf.multiply(t,t))))

y = bf.multiply(bf.multiply(2,t), bf.reciprocal(bf.subtract(1, bf.multiply(t,t))))



print(f"t: {t}")
print(f"-->  ({x}, {y})")

print(f"verify: 1 = {bf.subtract(bf.multiply(x,x),bf.multiply(y,y))}")

print(f"{bf.reciprocal(1)}\n\n\n")"""


bf = fieldmath.Zp(11)
d = 4
print(f"Must be true: {bf.multiply(d,d)} = {bf.add(1, bf.multiply(4,2))}")
e = bf.negate(1)

a = bf.multiply(e, d)
c = bf.multiply(e, bf.multiply(2, d))


n_one = bf.negate(1)

G = fieldmath.create_matrix([[a,1,1,1,1,1], [1,a,n_one,1,1,n_one],[1,n_one,a,n_one,1,1],[1,1,n_one,a,n_one,1],[1,1,1,n_one,a,n_one],[1,n_one,1,1,n_one,a]], bf)

print("G=\n", G, "\n\n\n")
print("kernel of G\n", G.kernel_space(), "\n\n\n")
print("G*G=\n", G*G, "\n\n\n")
print("c*G-G^2=\n",G*G-fieldmath.identity_n(6, bf, c)*G, "\n\n\n")