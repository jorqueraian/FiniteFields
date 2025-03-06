import networkx as nx
import os, sys
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import math
import itertools

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


def select_n_compat_vecs(all_vecs, n, search_space=None, vecs_chosen=[], b=1, lazy=False):
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
                if lazy:
                    return sclr_prod
                if isinstance(vec1.f, fieldmath.FieldWithInvolution):
                    return vec1.f.modulus_squared(sclr_prod)
                else:
                    ## Assume symmetric
                    return vec1.f.multiply(sclr_prod,sclr_prod)

            for i, add_vec in enumerate(search_space):
                #if i % (len(search_space)//20+1) != 0:
                #    continue
                compat_vecs = [j for ind, j in enumerate(search_space) if ind> i and scalar_product_mag(add_vec, j) == b ]
                yield from select_n_compat_vecs(all_vecs, n-1, compat_vecs, vecs_chosen+[add_vec])



def find_frame_pls(p, l, f, a, d, s, b=1, all_vecs=None, scalar_product_gram_discr=None):
    keys = ["p", "l", "f", "s", "s^2", "a", "d", "n"]
    #for (p, l, f, s, ss, a, d, n) in iter_numbers():
    if all_vecs is None:
        all_vecs = mat_of_non_zero_vectors(f, d, a, scalar_product_gram_discr)

    def scalar_product_mag(v1, v2, lazy=False):
        vec1 = all_vecs.get_sub_matrix_from_cols([v1])
        vec2 = all_vecs.get_sub_matrix_from_cols([v2])
        sclr_prod = (vec1.adjoint()*vec2).get(0,0)
        if lazy:
            return sclr_prod
        if isinstance(vec1.f, fieldmath.FieldWithInvolution):
            return vec1.f.modulus_squared(sclr_prod)
        else:
            ## Assume symmetric
            return vec1.f.multiply(sclr_prod,sclr_prod)

    for init_2_vecs in select_n_compat_vecs(all_vecs, 2, b=b, lazy=True):
        if init_2_vecs is None:
            continue
        else:
            # work from here
            verts = [j for j in range(all_vecs.columns) if j not in init_2_vecs and scalar_product_mag(j, j, True) == a and scalar_product_mag(j, init_2_vecs[0], lazy=True) == b and scalar_product_mag(j, init_2_vecs[1]) == b ]
            
            G = nx.Graph()
            G.add_nodes_from(verts)

            for u in verts:
                for v in verts:
                    if u != v and scalar_product_mag(u, v) == b:
                        G.add_edge(u,v)

            for K in nx.enumerate_all_cliques(G):
                if len(K) <= d-2:
                    continue
                maybe_frame = all_vecs.get_sub_matrix_from_cols(init_2_vecs+K)
                #if maybe_frame.rank() != d:
                #    continue

                gram = (maybe_frame.adjoint()*maybe_frame)

                if maybe_frame.rank() != gram.rank():
                    continue
                        
                # Now we know its a frame

                (is_ab, (true_a, b)) = frame_thy.is_equiangular(gram_mat=gram)
                #if not is_ab or b != 1:
                #    continue
                (is_c, c) = frame_thy.is_tight(gram_mat=gram)
                if not is_c:
                    continue
                # Now we know its (a,1,c)-ETF
                print(".", end="", flush=True)
                # print("ETF found\n", maybe_frame, maybe_frame_ind)
                # Now we want to hint for simplices

                binder = frame_thy.contains_simplex(s,maybe_frame)
                if len(binder) != 0:
                    return maybe_frame, binder
    return False


find_frame_pls(3, None, fieldmath.Zp(3), 0,7, 3, scalar_product_gram_discr=1)
print("2")
find_frame_pls(3, None, fieldmath.Zp(3), 0,7, 3, scalar_product_gram_discr=2)
print("3")
find_frame_pls(3, None, fieldmath.Zp(3), 0,9, 3, scalar_product_gram_discr=2)
print("4")
find_frame_pls(3, None, fieldmath.Zp(3), 0,9, 3, scalar_product_gram_discr=1)