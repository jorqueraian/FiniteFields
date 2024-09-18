# For field implementation look at fieldmath.py
import fieldmath


f11 = fieldmath.Zp(11)

print(f11.multiply(4,4), " = 5? You have found the square root of 5")
phi=f11.multiply(f11.add(1,4), f11.reciprocal(2))


G = fieldmath.create_matrix(
    [
        [0,    0,                1,    1,                phi,  f11.negate(phi)],
        [1,    1,                phi,  f11.negate(phi),  0,    0,],
        [phi,  f11.negate(phi),  0,    0,                1,    1,],
    ], f11)


print("G=\n", G, "\n\n\n")
print("|G^t*G|^2=\n", (G.conj_transpose()*G).modulus_squared_of_entries(), "\n\n\n")
print("G*G^t=\n", G*G.conj_transpose(), "\n\n\n")


"""columns_to_use = [0,1,2]
G_sub = G.get_sub_matrix_from_cols(columns_to_use)
print("G sub=\n", G_sub, "\n\n\n")
print("kernel of G sub\n", G_sub.kernel_space(), "\n\n\n")
print("|G_sub^t*G_sub|^2=\n", (G_sub.conj_transpose()*G_sub).modulus_squared_of_entries(), "\n\n\n")
print("G_sub*G_sub^t*G_sub*2=\n", G_sub*G_sub.conj_transpose()*G_sub*2, "\n\n\n")
print("In this case G_sub is tight with c'=3\n\n\n")"""


"""
columns_to_use2 = [3,4,5,6,7,8]
G_sub_comp = G.get_sub_matrix_from_cols(columns_to_use2)
print("G sub comp=\n", G_sub_comp, "\n\n\n")
print("kernel of G sub\n", G_sub_comp.kernel_space(), "\n\n\n")
print("|G_sub_comp^t*G_sub_comp|^2=\n", (G_sub_comp.conj_transpose()*G_sub_comp).modulus_squared_of_entries(), "\n\n\n")
#print("G_sub_comp^t*G_sub_comp=\n", G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
#print("(G_sub_comp^t*G_sub_comp)^2=\n", G_sub_comp.conj_transpose()*G_sub_comp*G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
print("G_sub_comp*G_sub_comp^t*G_sub_comp=\n", G_sub_comp*G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
print("In this case G_sub_comp is tight with c'=3")
"""