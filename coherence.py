import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import math
import itertools

def mat_of_non_zero_vectors(fld, d):
    num_fld_elms = fld.size
    num_vecs = num_fld_elms**d

    assert num_vecs < 10000, "You are joking right? do you want your computer to explode?"

    big_mat = [[] for _ in range(num_vecs)]

    fld_elms = [x for x in fld.iter_elems()]
    for i, vec in enumerate(itertools.product(fld_elms, repeat=d)):
        big_mat[i] = list(vec)
    
    return fieldmath.create_matrix(big_mat, fld).transpose()

def coherence_sqrd_mat(mat):
    f = mat.f
    gram = (mat.conj_transpose()*mat)
    diag_gram = fieldmath.create_matrix([[gram.get(i,i) for i in range(gram.columns)]], mat.f)
    coherence_mat = gram.modulus_squared_of_entries()
    
    for r in range(coherence_mat.rows):
        for c in range(coherence_mat.columns):
            grr = gram.get(r,r)
            gcc = gram.get(c,c)

            if f.equals(grr, f.zero()) or f.equals(gcc, f.zero()):
                coherence_mat.set(r, c, None)
            else:
                coherence_mat.set(r, c, f.divide(coherence_mat.get(r,c), f.multiply(grr,gcc)))
    return coherence_mat, diag_gram

def find_equiangular_systems_from_coherence_mat_3vecs(coherence_mat, vec_mags, all_vecs, coherence_val=None):
    f = coherence_mat.f
    sqrs = list(set([f.print_elm(f.multiply(x, x)) for x in f.iter_elems()]))
    def is_sqr(x):
        if f.equals(x, f.zero()):
            return False
        if (f.print_elm(x) in sqrs):
            return True
        else:
            return False
    
    equiang_syss = []
    special_equiang_syss = []
    for j1 in range(coherence_mat.columns):
        for j2 in range(j1+1, coherence_mat.columns):
            for j3 in range(j2+1, coherence_mat.columns):
                if coherence_mat.get(j1,j2) is not None and coherence_mat.get(j1,j3) is not None and coherence_mat.get(j2,j3) is not None:
                    if f.equals(coherence_mat.get(j1,j2),coherence_mat.get(j2,j3)) and f.equals(coherence_mat.get(j1,j3),coherence_mat.get(j2,j3)): # and not f.equals(coherence_mat.get(j1,j3),f.zero()):
                        if all_vecs.get_sub_matrix_from_lists([i for i in range(all_vecs.rows)], [j1,j2,j3]).rank() != 1:
                            equiang_syss.append((j1,j2,j3))
                            if not f.equals(vec_mags.get(0, j1),vec_mags.get(0, j2)) or not f.equals(vec_mags.get(0, j1),vec_mags.get(0, j3)):
                                if not is_sqr(vec_mags.get(0, j1)) or not is_sqr(vec_mags.get(0, j2)) or not is_sqr(vec_mags.get(0, j3)):
                                    if is_sqr(vec_mags.get(0, j1)) or is_sqr(vec_mags.get(0, j2)) or is_sqr(vec_mags.get(0, j3)):
                                        special_equiang_syss.append((j1,j2,j3))
    return equiang_syss, special_equiang_syss


def find_equiangular_frames_from_coherence_mat(coherence_mat, vec_mags, all_vecs, num_vecs, coherence_val=None):
    f = coherence_mat.f
    sqrs = list(set([f.print_elm(f.multiply(x, x)) for x in f.iter_elems()]))
    def is_sqr(x):
        if f.equals(x, f.zero()):
            return False
        if (f.print_elm(x) in sqrs):
            return True
        else:
            return False
    
    equiang_syss = []
    special_equiang_syss = []
    for j1 in range(coherence_mat.columns):
        for j2 in range(j1+1, coherence_mat.columns):
            for j3 in range(j2+1, coherence_mat.columns):
                if coherence_mat.get(j1,j2) is not None and coherence_mat.get(j1,j3) is not None and coherence_mat.get(j2,j3) is not None:
                    if f.equals(coherence_mat.get(j1,j2),coherence_mat.get(j2,j3)) and f.equals(coherence_mat.get(j1,j3),coherence_mat.get(j2,j3)):
                        if all_vecs.get_sub_matrix_from_lists([i for i in range(all_vecs.rows)], [j1,j2,j3]).rank() != 1:
                            equiang_syss.append((j1,j2,j3))
                            if not f.equals(vec_mags.get(0, j1),vec_mags.get(0, j2)) or not f.equals(vec_mags.get(0, j1),vec_mags.get(0, j3)):
                                if not is_sqr(vec_mags.get(0, j1)) or not is_sqr(vec_mags.get(0, j2)) or not is_sqr(vec_mags.get(0, j3)):
                                    if is_sqr(vec_mags.get(0, j1)) or is_sqr(vec_mags.get(0, j2)) or is_sqr(vec_mags.get(0, j3)):
                                        special_equiang_syss.append((j1,j2,j3))
    return equiang_syss, special_equiang_syss

            
#p=3# also try 47, 59, 61. # already tried 37, 47 and lower primes
#fld = fieldmath.Zp(p) #fieldmath.FieldExtension(fieldmath.Zp(p), [p-3,0,1], 1)
#all_vecs = mat_of_non_zero_vectors(fld, 3)
#coherencesqrd_mat, vec_mags = coherence_sqrd_mat(all_vecs)
#l1, l2 = find_equiangular_systems_from_coherence_mat_3vecs(coherencesqrd_mat, vec_mags, all_vecs)
#print(l2)


