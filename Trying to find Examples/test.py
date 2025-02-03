# For field implementation look at fieldmath.py
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)
d = "-1"

#print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

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


Phi= hessa_sic
G=(Phi.conj_transpose())*Phi
is_G_frm, (rnk_G, disc_G) = frame_thy.is_frame(None, G, False)
isetf, (a,b,c_G) = frame_thy.is_etf(None, G, is_G_frm)

print(f"Phi is an ({a},{b},{c_G})-ETF of a {rnk_G}-dimensional non-degenerate space")

subframe_vecs = [0,1,3]
Phi_simplex = Phi.get_sub_matrix_from_cols(subframe_vecs)
G_simplex=(Phi_simplex.conj_transpose())*Phi_simplex
is_simp_frm, (rnk_sim, disc_simp) = frame_thy.is_frame(None, G_simplex, False)
if rnk_sim == Phi_simplex.rank():
    isimpetf, (a,b,c_simp) = frame_thy.is_etf(None, G_simplex, is_simp_frm)
    if isimpetf:
        print(f"The vectors {subframe_vecs} form a sub-({a},{b},{c_simp})-ETF for its {rnk_sim}-dimension span")
    else: 
        print(f"The vectors {subframe_vecs} span a non-degerate space but do not form a sub-ETF for its {rnk_sim}-dimension span")
else:
    print(f"The vectors {subframe_vecs} do not form a frame and span a degenerate space.")

"""
print("Phi=\n", Phi, "\n\n\n")
print("Phi^t*Phi=\n", (Phi.conj_transpose())*Phi, "\n\n\n")
print("Phi*Phi^t=\n", Phi*(Phi.conj_transpose()), "\n\n\n")


print(f25.multiply("1+1a",f25.multiply("1+1a","1+1a")))"""
