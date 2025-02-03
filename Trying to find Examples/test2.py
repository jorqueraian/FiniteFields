import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math

f125 = fieldmath.FieldExtension(fieldmath.Zp(5),[1,1,0,1])
f25 = fieldmath.FieldExtension(fieldmath.Zp(5),[3,0,1])



p=7
d=13
n=26
s=5
conf_fld = f25

fld = fieldmath.Zp(p)
G = frame_thy.create_gram_of_d_2d_etf_from_conference_mat(
    frame_thy.create_conference_matrix(
        conf_fld, 
        fieldmath.Zp(p)
    ), 
    fld, 
    s, 
    True
    )

compute_discrim = False

is_G_frm, (rnk_G, disc_G) = frame_thy.is_frame(None, G, compute_discrim)
isetf, (a,b,c_G) = frame_thy.is_etf(None, G, is_G_frm)

assert is_G_frm and isetf, "yikes"
# We can test some possible frame vectors
its = 0
print("There are: ", (math.comb(G.rows, s+1)), "things to check")
# Note I do believe $s=13 as well
for sub_frame_G, sub_frame_cols in fieldmath.iter_sub_mats(G, s+1, randomize=False):
    sub_frame_G_cols = G.get_sub_matrix_from_cols(sub_frame_cols)
    if its % 1000 == 999:
        print(".", end='', flush=True)
    its+=1
    if sub_frame_G_cols.rank() != s:
        continue
    
    is_Gsub_frm, (rnk_Gsub, disc_Gsub) = frame_thy.is_frame(None, sub_frame_G, compute_discrim)
    is_tgt, c = frame_thy.is_tight(None, sub_frame_G)
    if is_tgt and is_Gsub_frm and sub_frame_G_cols.rank() == rnk_Gsub:
        print(f"\nYo! the frame vectors {list(sub_frame_cols)} is a ({a},{b},{c})-ETF, therefore a regular {s}-simplex!")
        break
print("\nSad stuff, no regular simplex here, moving on to the next")