# For field implementation look at fieldmath.py
import fieldmath


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)
d = "-1"

print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

hessa_sic = fieldmath.create_matrix(
    [
        ["1","1","1",         "-1","-1","-1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",     "-1","-1a","-1a^2"],
        ["-1","-1a^2","-1a",  "0","0","0",         "1","1a^2","1a"   ]
    ], f25)*"1"

emilys_fav_frame = fieldmath.create_matrix(
    [
        ["1","1","1",         "1","1","1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",  "1","1a","1a^2"   ],
        ["1","1a^2","1a",     "0","0","0",      "1","1a^2","1a"   ]
    ], f25)

G_other = fieldmath.create_matrix(
    [
        [1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,-1,-1,-1,-1,-1,-1],
        [1,-1,-1,-1,1,1,1,-1,-1,-1],
        [-1,1,-1,-1,1,-1,-1,1,1,-1],
        [-1,-1,1,-1,-1,1,-1,1,-1,1],
        [-1,-1,-1,1,-1,-1,1,-1,1,1],
    ], fieldmath.Zp(3))


G= hessa_sic
print("G=\n", G, "\n\n\n")
print("|G^t*G|^2=\n", (G.conj_transpose()*G).modulus_squared_of_entries(), "\n\n\n")
print("G*G^t=\n", G*G.conj_transpose(), "\n\n\n")

columns_to_use = [0,1,2]
G_sub = G.get_sub_matrix_from_cols(columns_to_use)
print("G sub=\n", G_sub, "\n\n\n")
print("kernel of G sub\n", G_sub.kernel_space(), "\n\n\n")
print("|G_sub^t*G_sub|^2=\n", (G_sub.conj_transpose()*G_sub).modulus_squared_of_entries(), "\n\n\n")
print("G_sub*G_sub^t*G_sub*2=\n", G_sub*G_sub.conj_transpose()*G_sub*2, "\n\n\n")
print("In this case G_sub is tight with c'=3\n\n\n")


columns_to_use2 = [3,4,5,6,7,8]
G_sub_comp = G.get_sub_matrix_from_cols(columns_to_use2)
print("G sub comp=\n", G_sub_comp, "\n\n\n")
print("kernel of G sub\n", G_sub_comp.kernel_space(), "\n\n\n")
print("|G_sub_comp^t*G_sub_comp|^2=\n", (G_sub_comp.conj_transpose()*G_sub_comp).modulus_squared_of_entries(), "\n\n\n")
#print("G_sub_comp^t*G_sub_comp=\n", G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
#print("(G_sub_comp^t*G_sub_comp)^2=\n", G_sub_comp.conj_transpose()*G_sub_comp*G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
print("G_sub_comp*G_sub_comp^t*G_sub_comp=\n", G_sub_comp*G_sub_comp.conj_transpose()*G_sub_comp, "\n\n\n")
print("In this case G_sub_comp is tight with c'=3")
