# For field implementation look at fieldmath.py
import fieldmath


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)
d = "1+1a"

print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

G = fieldmath.create_matrix(
    [
        ["1","1","1",      "0","0","0",      "1","1","1"],
        ["0","0","0",      "1","1a","1a^2",  "1","1a","1a^2"],
        ["1","1a^2","1a",  "1","1a^2","1a",  "0","0","0"]
    ], f25)

G1 = fieldmath.create_matrix(
    [
        ["1","1","1"],
        ["1","1","1"],
        ["1","1a^2","1a"]
    ], f25)

print("G=\n", G, "\n\n\n")
print("kernel of G\n", G.kernel_space(), "\n\n\n")
print("G=\n", G, "\n\n\n")
print("G^t*G=\n", G.conj_transpose()*G, "\n\n\n")
print("G*G^t=\n", G*G.conj_transpose(), "\n\n\n")