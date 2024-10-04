# For field implementation look at fieldmath.py
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath


f7 = fieldmath.Zp(7)
f11 = fieldmath.Zp(11)
#print(f11.reciprocal(4))

# now we want to bake in sqrt(5), as it doesnt exist
f49 = fieldmath.FieldExtension(f7, [-5,0,1], involution_pow=None)
#print(f49.multiply("1a", "3a"))

conference_49 = fieldmath.create_matrix(
    [
        [0,1,1,1,1,1  ],
        [1,0,1,-1,-1,1],
        [1,1,0,1,-1,-1],
        [1,-1,1,0,1,-1],
        [1,-1,-1,1,0,1],
        [1,1,-1,-1,1,0]
    ], f49)
conference_11 = fieldmath.create_matrix(
    [
        [0,1,1,1,1,1  ],
        [1,0,1,-1,-1,1],
        [1,1,0,1,-1,-1],
        [1,-1,1,0,1,-1],
        [1,-1,-1,1,0,1],
        [1,1,-1,-1,1,0]
    ], f11)

G_49=conference_49*"3a" + fieldmath.identity_n(6, f49)
G_11=conference_11*3 + fieldmath.identity_n(6, f11)


print("over F_49 G=\n", G_49, "\n\n\n")
print("over F_49 G*G-2*G=\n", G_49*G_49-G_49*"2", "\n\n\n")
print("over F_49 G is the gram matrix of a tight frame with c=2\n\n\n")

#print("over F_11 G=\n", G_11, "\n\n\n")
#print("over F_11 G*G-2*G=\n", G_11*G_11-G_11*2, "\n\n\n")
#print("over F_11 G is the gram matrix of a tight frame with c=2\n\n\n")


# Now we want to find a a basic sub matrix, ie a basis for the span of G

"""
columns_to_use = [0,1,2]
print("kernel of G\n", G_49.kernel_space(), "\n\n\n")
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