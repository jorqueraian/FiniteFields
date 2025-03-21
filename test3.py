import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import math
import itertools
from coherence import mat_of_non_zero_vectors, coherence_sqrd_mat

def find_ETFs_from_coherence_mat(coherence_mat, vec_mags, all_vecs, num_vecs, coherence_val=None, s=None):
    f = coherence_mat.f
    sqrs = list(set([f.print_elm(f.multiply(x, x)) for x in f.iter_elems()]))
    def is_sqr(x):
        if f.equals(x, f.zero()):
            return False
        if (f.print_elm(x) in sqrs):
            return True
        else:
            return False
    
    #print(f"there are {math.comb(coherence_mat.rows, num_vecs)} combinations to check. That probably way too many but ill give it a go!")
    for first_vec in range(0, coherence_mat.rows):
        fv_norm = vec_mags.get(0,first_vec)
        if ((s % f.p) == 0 and f.equals(fv_norm, f.zero())) or ((s % f.p) != 0 and not f.equals(fv_norm, f.zero())):
            compatible_vecs_inds = [i for i in range(first_vec, coherence_mat.rows) if i != first_vec and vec_mags.get(0,i) == fv_norm]

            print(f"there are {math.comb(len(compatible_vecs_inds), num_vecs-1)} combinations to check, for subsets containing first vector. ({first_vec+1}/{coherence_mat.rows})")
            for parsubset in itertools.combinations(compatible_vecs_inds, num_vecs-1):
                # Check if ETF
                subset = [first_vec] + list(parsubset)
                # norms = vec_mags.get_sub_matrix_from_cols(subset)
                if True: # all([f.equals(norms.get(0,0), norms.get(0,i)) for i in range(1, norms.columns)]):
                    phi = all_vecs.get_sub_matrix_from_cols(subset)
                    gram_mat = (phi.conj_transpose()*phi)
                    
                    (is_ab_equi,(a,b)) = frame_thy.is_equiangular(gram_mat=gram_mat)
                    if is_ab_equi:
                        #if not f.equals(f.multiply(a,a), f.multiply(f.multiply(s,s),b)):
                        #    pass
                        #print("We got an Equiangular System")
                        if frame_thy.is_frame(gram_mat=gram_mat)[0]:
                            #print("We got a frame!")
                            if frame_thy.is_tight(gram_mat=gram_mat)[0]:
                                #print("We got an ETF!")
                                spk = frame_thy.spark(phi)
                                if spk > 2:
                                    print(spk)
                                    print(subset)
                                #yield phi, gram_mat

# We want to find an example here with a^2!=b is orthog geom, spark=3. They dont exist
# what about spark 4?
# with s=2 internal simplex
n=6
d=3
p=7 # also try 7
s=3 #try 2,
fld = fieldmath.Zp(p) #fieldmath.FieldExtension(fieldmath.Zp(p), [p-3,0,1], 1)
all_vecs = mat_of_non_zero_vectors(fld, 3)
coherencesqrd_mat, vec_mags = coherence_sqrd_mat(all_vecs)
find_ETFs_from_coherence_mat(coherencesqrd_mat, vec_mags, all_vecs, num_vecs=n, s=s)
print("DOne!")