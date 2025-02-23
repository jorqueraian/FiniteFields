import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math
import itertools


def iter_numbers():
    keys = ["p", "l", "f", "s", "s^2", "a", "d", "n"]
    results = []
    for p in [x for x in range(5,30) if ians_numers.isprime(x)]:
        if p % 4 != 1:
            continue

        f = fieldmath.Zp(p)
        sqrs = [(x,f.print_elm(f.multiply(x, x))) for x in f.iter_elems()]
        def square_root(x):
            return next((srx[0] for srx in sqrs if srx[1]==x), None)
        
        for d in range(5, 12):
            for n in range(d+2, d*d):
                if (n-d) % p == 0:
                    continue
                welch_ss = f.divide(f.multiply(d, n-1) , n-d)
                s = square_root(welch_ss)
                
                if s is None or s % p == 0 or s+1 % p == 0:
                    continue
                
                # s needs to be the negative of a
                a = f.negate(s)
                if s != 1 and s < d:
                    yield (p, None, f, s, welch_ss, a, d, n)
                if a != 1 and a < d:
                    yield (p, None, f, a, welch_ss, s, d, n)


def mat_of_non_zero_vectors(fld, d, mag):
    num_fld_elms = fld.size
    max_num_vecs = num_fld_elms**d

    assert max_num_vecs/num_fld_elms < 1000000000, "You are joking right? do you want your computer to explode?"

    big_mat = []# [[] for _ in range(num_vecs)]

    fld_elms = [x for x in fld.iter_elems()]
    for i, vec in enumerate(itertools.product(fld_elms, repeat=d)):
        mag_sqrd = sum([v**2 for v in vec])
        if mag_sqrd == mag:
            big_mat.append(list(vec))
    
    return fieldmath.create_matrix(big_mat, fld).transpose()


def select_n_compat_vecs(all_vecs, all_vecs_gram, n, search_space=None, vecs_chosen=[]):
    if n == 0:
        yield vecs_chosen
    
    if search_space is None:
        search_space = [j for j in range(all_vecs.columns)]

    if len(search_space) < n:
        yield None

    f = all_vecs.f

    for i, add_vec in enumerate(search_space):
        compat_vecs = [j for j in search_space if j> i and f.multiply(all_vecs_gram.get(add_vec, j), all_vecs_gram.get(j, add_vec)) == 1 ]
        yield from select_n_compat_vecs(all_vecs, all_vecs_gram, n-1, compat_vecs, vecs_chosen+[add_vec])


def check_thy_numbers(p, l, f,s, ss, a, d, n):
    keys = ["p", "l", "f", "s", "s^2", "a", "d", "n"]
    #for (p, l, f, s, ss, a, d, n) in iter_numbers():
    all_vecs = mat_of_non_zero_vectors(f, d, a)
    all_vecs_gram = (all_vecs.conj_transpose()*all_vecs)
    for maybe_frame_ind in select_n_compat_vecs(all_vecs, all_vecs_gram, n):
        if maybe_frame_ind is None:
            continue
        else:
            maybe_frame = all_vecs.get_sub_matrix_from_cols(maybe_frame_ind)
            gram = (maybe_frame.conj_transpose()*maybe_frame)

            if maybe_frame.rank() != d:
                continue
                    
            # Now we know its a frame

            (is_ab, (true_a, b)) = frame_thy.is_equiangular(gram_mat=gram)
            if not is_ab or b != 1:
                continue
            (is_c, c) = frame_thy.is_tight(gram_mat=gram)
            if not is_c:
                continue
            # Now we know its (a,1,c)-ETF
            # Now we want to hint for simplices

            for maybe_simplex_ind in itertools.combinations(range(maybe_frame.columns), s+1):
                maybe_simplex = maybe_frame.get_sub_matrix_from_cols(maybe_simplex_ind)
                simpl_gram = (maybe_simplex.conj_transpose()*maybe_simplex)

                if maybe_simplex.rank() != s and simpl_gram.rank() != s:
                    continue
                
                # Now we know its a frame

                (is_c_simp, c_simp) = frame_thy.is_tight(gram_mat=simpl_gram)
                    
                if not is_c_simp:
                    continue

                return maybe_frame, maybe_simplex, maybe_frame_ind, maybe_simplex_ind
    return False


#for row in iter_numbers():
#    print(row)

# NO GOs
#(5, None, <fieldmath.Zp object at 0x000002671875DBE0>, 3, 4, 2, 8, 59)
# (5, None, <fieldmath.Zp object at 0x000002671875DBE0>, 2, 4, 3, 7, 13)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 7, 15)
#(5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 3, 4, 2, 7, 18)
#(5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 7, 20)

## TRY THESE
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 8, 10)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 2, 4, 3, 8, 14)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 2, 4, 3, 8, 19)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 9, 15)

# (13, None, <fieldmath.Zp object at 0x000001B4F5776D50>, 4, 3, 9, 5, 8)
# (13, None, <fieldmath.Zp object at 0x00000267187A6D50>, 3, 9, 10, 5, 10)
# (13, None, <fieldmath.Zp object at 0x000001B4F5776D50>, 5, 12, 8, 6, 11)

for (p, _, f, s, welch_ss, a, d, n) in iter_numbers():
    print(p, _, f, s, welch_ss, a, d, n)
    print(check_thy_numbers(p, None, f, s, welch_ss, a, d, n))

