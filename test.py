# For field implementation look at fieldmath.py
import fieldmath


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)
d = "-1"

print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

G = fieldmath.create_matrix(
    [
        ["1","1","1",         "-1","-1","-1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",     "-1","-1a","-1a^2"],
        ["-1","-1a^2","-1a",  "0","0","0",         "1","1a^2","1a"   ]
    ], f25)*"1"

G1 = fieldmath.create_matrix(
    [
        ["1","1","1",         "1","1","1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",     "1","1a","1a^2"],
        ["1","1a^2","1a",  "0","0","0",         "1","1a^2","1a"   ]
    ], f25)

G2 = fieldmath.create_matrix(
    [
        ["1","1","1"],
        ["1","1","1"],
        ["1","1a^2","1a"]
    ], f25)

print("G=\n", G, "\n\n\n")
columns_to_use = [1,4,7]
print("G sub=\n", G.get_sub_matrix_from_cols(columns_to_use), "\n\n\n")
print("kernel of G sub\n", G.get_sub_matrix_from_cols(columns_to_use).kernel_space(), "\n\n\n")

#print("G=\n", G, "\n\n\n")
#print("|G^t*G|^2=\n", (G.conj_transpose()*G).modulus_squared_of_entries(), "\n\n\n")
#print("G*G^t=\n", G*G.conj_transpose(), "\n\n\n")