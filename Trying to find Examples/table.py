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


def mat_of_non_zero_vectors(fld, d, vec_mag_sqrd, output_gram_discr=None):
    num_fld_elms = fld.size
    max_num_vecs = num_fld_elms**d

    if output_gram_discr is None:
        output_gram = [1]*d
    else:
        output_gram = [1]*(d-1) + [output_gram_discr]

    if max_num_vecs/num_fld_elms > 100000000:  #, "You are joking right? do you want your computer to explode?"
        print("You are joking right? do you want your computer to explode? Im going to just pretend you were joking and ignore that")
        return None

    big_mat = []# [[] for _ in range(num_vecs)]

    fld_elms = [x for x in fld.iter_elems()]
    for i, vec in enumerate(itertools.product(fld_elms, repeat=d)):
        mag_sqrd = fieldmath.field_sum(fld, [fld.multiply(fld.multiply(v,v), og) for v, og in zip(vec, output_gram)])
        if mag_sqrd == vec_mag_sqrd:
            big_mat.append(list(vec))
    
    return fieldmath.create_matrix(big_mat, fld, row_space_gram=output_gram).transpose()


def select_n_compat_vecs(all_vecs, n, search_space=None, vecs_chosen=[], b=1):
    if n == 0:
        yield vecs_chosen
    else:
        if search_space is None:
            search_space = [j for j in range(all_vecs.columns)]

        if len(search_space) < n:
            yield None
        else:

            f = all_vecs.f

            def scalar_product_mag(v1, v2):
                vec1 = all_vecs.get_sub_matrix_from_cols([v1])
                vec2 = all_vecs.get_sub_matrix_from_cols([v2])
                sclr_prod = (vec1.adjoint()*vec2).get(0,0)
                if isinstance(vec1.f, fieldmath.FieldWithInvolution):
                    return vec1.f.modulus_squared(sclr_prod)
                else:
                    ## Assume symmetric
                    return vec1.f.multiply(sclr_prod,sclr_prod)

            for i, add_vec in enumerate(search_space):
                compat_vecs = [j for ind, j in enumerate(search_space) if ind> i and scalar_product_mag(add_vec, j) == b ]
                yield from select_n_compat_vecs(all_vecs, n-1, compat_vecs, vecs_chosen+[add_vec])


def check_thy_numbers(p, l, f, s, ss, a, d, n, b=1, all_vecs=None, scalar_product_gram_discr=None):
    keys = ["p", "l", "f", "s", "s^2", "a", "d", "n"]
    #for (p, l, f, s, ss, a, d, n) in iter_numbers():
    if all_vecs is None:
        all_vecs = mat_of_non_zero_vectors(f, d, a, scalar_product_gram_discr)

    for maybe_frame_ind in select_n_compat_vecs(all_vecs, n, b=b):
        if maybe_frame_ind is None:
            continue
        else:
            maybe_frame = all_vecs.get_sub_matrix_from_cols(maybe_frame_ind)
            if maybe_frame.rank() != d:
                continue

            gram = (maybe_frame.adjoint()*maybe_frame)

                    
            # Now we know its a frame

            #(is_ab, (true_a, b)) = frame_thy.is_equiangular(gram_mat=gram)
            #if not is_ab or b != 1:
            #    continue
            (is_c, c) = frame_thy.is_tight(gram_mat=gram)
            if not is_c:
                continue
            # Now we know its (a,1,c)-ETF
            print(".", end="", flush=True)
            # print("ETF found\n", maybe_frame, maybe_frame_ind)
            # Now we want to hint for simplices

            for maybe_simplex_ind in itertools.combinations(range(maybe_frame.columns), s+1):
                maybe_simplex = maybe_frame.get_sub_matrix_from_cols(maybe_simplex_ind)
                simpl_gram = (maybe_simplex.adjoint()*maybe_simplex)

                if maybe_simplex.rank() != s and simpl_gram.rank() != s:
                    continue
                
                # Now we know its a frame

                (is_c_simp, c_simp) = frame_thy.is_tight(gram_mat=simpl_gram)
                    
                if not is_c_simp:
                    continue

                return maybe_frame, maybe_simplex, maybe_frame_ind, maybe_simplex_ind
    return False


def check_thy_numbers_loop(initial_p=0, initial_d=0, initial_n=0):
    old_p = 0
    old_d = 0
    max_n = 0
    old_a = 0
    bad_as = []
    all_vecs = None
    print("p, s, s^2, a, d, n")
    for (p, _, f, s, welch_ss, a, d, n) in iter_numbers():
        if p < initial_p or d < initial_d or n < initial_n:
            continue

        print(p, s, welch_ss, a, d, n,)
        #print(check_thy_numbers(p, None, f, s, welch_ss, a, d, n))
        if all_vecs is None or old_d != d or old_p != p or old_a != a:
            all_vecs = mat_of_non_zero_vectors(f, d, a)
            if all_vecs is None:
                continue

            max_n = d**2
            if old_d != d or old_p != p:
                bad_as = []
        old_p = p
        old_d = d
        
        # Changes based on n
        # But can n vectors even be found?
        if n > max_n and a in bad_as:
            continue
        
        vecs_found = False
        for maybe_frame_ind in select_n_compat_vecs(all_vecs, n):
            if maybe_frame_ind is None:
                continue
            else:
                vecs_found = True
                maybe_frame = all_vecs.get_sub_matrix_from_cols(maybe_frame_ind)
                if maybe_frame.rank() != d:
                    continue

                gram = (maybe_frame.adjoint()*maybe_frame)
                        
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
                    simpl_gram = (maybe_simplex.adjoint()*maybe_simplex)

                    if maybe_simplex.rank() != s and simpl_gram.rank() != s:
                        continue
                    
                    # Now we know its a frame

                    (is_c_simp, c_simp) = frame_thy.is_tight(gram_mat=simpl_gram)
                        
                    if not is_c_simp:
                        continue

                    return maybe_frame, maybe_simplex, maybe_frame_ind, maybe_simplex_ind
        if not vecs_found:
            max_n = n
            bad_as.append(a)


#for row in iter_numbers():
#    print(row)


#(5, None, <fieldmath.Zp object at 0x000002671875DBE0>, 3, 4, 2, 8, 59)
# (5, None, <fieldmath.Zp object at 0x000002671875DBE0>, 2, 4, 3, 7, 13)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 7, 15)
#(5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 3, 4, 2, 7, 18)
#(5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 7, 20)

## TRY THESE

#check_thy_numbers(5, None, fieldmath.Zp(5),2, 4, 3, 7, 13, b=1, scalar_product_gram_discr=3)

# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 8, 10)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 2, 4, 3, 8, 14)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 2, 4, 3, 8, 19)
# (5, None, <fieldmath.Zp object at 0x000001B4F572DBE0>, 4, 1, 1, 9, 15)

# (13, None, <fieldmath.Zp object at 0x000001B4F5776D50>, 4, 3, 9, 5, 8)
# (13, None, <fieldmath.Zp object at 0x00000267187A6D50>, 3, 9, 10, 5, 10)
# (13, None, <fieldmath.Zp object at 0x000001B4F5776D50>, 5, 12, 8, 6, 11)

#for (p, _, f, s, welch_ss, a, d, n) in iter_numbers():
#    print(p, _, f, s, welch_ss, a, d, n)
#    print(check_thy_numbers(p, None, f, s, welch_ss, a, d, n))

#print(check_thy_numbers_loop(initial_p=5, initial_d=0))

#check_thy_numbers(5, None, fieldmath.Zp(5),3, 4, 2, 10, 45)

#Sanity

#all_vecs = mat_of_non_zero_vectors(fieldmath.Zp(5), 11, 2)
#all_vecs_gram = (all_vecs.adjoint()*all_vecs)
#for maybe_frame_ind in select_n_compat_vecs(all_vecs, all_vecs_gram, 56):
#    if maybe_frame_ind is not None:
#        print("OK verified")

#print("Ok bad things")




check_thy_numbers(3, None, fieldmath.Zp(3),3, 0, 0, 4, 10, b=1, scalar_product_gram_discr=2)
#print("DONE No simplex here! :()")

"""for p in [x for x in range(3,15) if ians_numers.isprime(x)]:

    f = fieldmath.Zp(p)
    sqrs = [(x,f.print_elm(f.multiply(x, x))) for x in f.iter_elems()]
    def square_root(x):
        return next((srx[0] for srx in sqrs if srx[1]==x), None)
        
    for a in f.iter_elems():
        for b in f.iter_elems():
            if f.equals(b,f.zero()):
                continue
            elif not f.equals(f.multiply(6,f.multiply(a,a)), f.multiply(4*9,b)):
                continue
            else:
                # dims line up
                check_thy_numbers(p, None, f, 0, 0, a, 4, 10, b=b)"""




#Here is the theory I want to build
# It is know that ETFs in the orthogonal geometry exists (In cases where real ETFs dont)
# however what is not know is if they exist in the real model, ie when the discr is a square.
# I want to be able to say that if a ETF in an orthogonal geometry exists then is must exists with some specific discriminant, or both.




### Here is a frame of 10 vecs over F_3^4
#  0  0  0  0  1  1  1  1  1  1  
#   0  0  1  1  0  0  1  1  2  2  
#   1  1  0  0  0  0  1  2  1  2  
#   1  2  1  2  1  2  0  0  0  0