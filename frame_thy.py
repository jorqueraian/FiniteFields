import fieldmath
import itertools
import math


def create_conference_matrix(construction_field, mat_field):
    # NOTE: THIS REQUIRES THAT THE USER KNOWS WHAT IMPORTS WOULD
    # SATISFY THE CONDITIONS OF THIS ALGORITHM
    # OR maybe add checks: finite field of p^k elements
    # will create a matrix that is p^k+1 X p^k+1, But it must be the 
    # case that $p^k$ is a sum of squares so either p not equiv 3 mod 4 or k is even.
    sqrs = list(set([construction_field.print_elm(construction_field.multiply(x, x)) for x in construction_field.iter_elems()]))
    
    def is_sqr(x):
        if construction_field.equals(x, construction_field.zero()):
            return mat_field.zero()
        if (construction_field.print_elm(x) in sqrs):
            return mat_field.one()
        else:
            return mat_field.negate(mat_field.one())

    C_lst = [[mat_field.zero()] + [mat_field.one()]*(construction_field.size)] + [
        [mat_field.one()]+
        [
            is_sqr(construction_field.subtract(y,x)) 
            for x in construction_field.iter_elems()
        ] 
        for y in construction_field.iter_elems()
    ]

    return fieldmath.create_matrix(C_lst, mat_field)


def create_gram_of_d_2d_etf_from_conference_mat(conference_mat, field, sqrt_2dminus1, promise_C_is_symmetric=False):
    assert conference_mat.f == field, "Sorry my code isnt smart enough to handle things like this"
    # this is an easy check but im lazy and the create_conf_mats function creates symmetric ones
    assert promise_C_is_symmetric, "Im sorry but you must promise me this one thing"

    return conference_mat*field.reciprocal(sqrt_2dminus1)+ fieldmath.identity_n(conference_mat.rows, field)


def is_equiangular(Phi=None, gram_mat=None, why_not=False):
    """Assumes equinorm, and then computes (a,b) parameters"""
    if Phi is None and gram_mat is None:
        assert False, "you need to specify an input"
    
    if gram_mat is None:
        gram_mat = (Phi.conj_transpose()*Phi)
    gram_mat_modulus_sqrd = gram_mat.modulus_squared_of_entries()
    a = gram_mat_modulus_sqrd.get(0, 0)  # Is this actually a^2? maybe
    b = gram_mat_modulus_sqrd.get(0, 1)

    # This is lazy but im ok with that.
    test_mat = fieldmath.identity_n(gram_mat.rows, gram_mat.f, gram_mat.f.subtract(a,b)) + fieldmath.Matrix(gram_mat.rows, gram_mat.rows, gram_mat.f, b)
    if (gram_mat_modulus_sqrd - test_mat).any():
        if why_not:
            print("Modulus Squared of Gram matrix:\n\n", gram_mat_modulus_sqrd)
        return False, (None, None)
    else:
        return True, (a,b)
    

def is_tight(Phi=None, gram_mat=None, why_not=False):
    if Phi is None and gram_mat is None:
        assert False, "you need to specify an input"
    
    if gram_mat is None:
        gram_mat = (Phi.conj_transpose()*Phi)

    gram_sqrd = gram_mat*gram_mat

    c = gram_mat.f.divide(gram_sqrd.get(0,0),gram_mat.get(0,0))

    if (gram_sqrd-(gram_mat*c)).any():
        if why_not:
            print("Its just not in the cards, sorry.")
        return False, None
    else:
        return True, c


def is_etf(Phi=None, gram_mat=None, do_you_promise_its_a_frame=True, why_not=False):
    
    assert do_you_promise_its_a_frame, "I need you to promise"

    is_equi, (a, b) = is_equiangular(Phi, gram_mat, why_not)
    if is_equi is False:
        return False, (None, None, None)

    is_tgt, c = is_tight(Phi, gram_mat, why_not)

    if is_tgt is False:
        return False, (a, b, None)
    else:
        return True, (a,b,c)


def is_frame(Phi=None, gram_mat=None, with_discr=True):
    # This is really hard as we need to check if the span in non degenerate,
    # And maybe something to do with the discriminant. But regardless, 
    # this comes down to computing some determinant, specifically of the gram matrix of the IP on the image of the frame
    # There are 3 things to check here 1) G
    assert gram_mat is not None, "Sorry i haven't done the other case yet"
    
    if gram_mat is not None:
        fld = gram_mat.f
    else:
        fld = Phi.f
    sqrs = list(set([fld.print_elm(fld.multiply(x, x)) for x in fld.iter_elems()]))
    
    def is_sqr(x):
        if fld.equals(x, fld.zero()):
            return fld.zero()
        if (fld.print_elm(x) in sqrs):
            return True
        else:
            return False
    
    if gram_mat is not None:
        if (gram_mat-gram_mat.conj_transpose()).any():
            return False, (None, None)
        rank = gram_mat.rank()

        if with_discr is False:
            return True, (rank, None)

    # only need to do this if the field automorphism is trivial
    discriminant = None
    total_columns_set = {k for k in range(gram_mat.columns)}

    for columns_lst in itertools.combinations(total_columns_set, rank):
        selected_columns = gram_mat.get_sub_matrix_from_cols(columns_lst)
        if selected_columns.rank() != rank:
            continue

        sub_mat = gram_mat.get_sub_matrix_from_lists(columns_lst, columns_lst)
        mat_det = sub_mat.det()
        discriminant = is_sqr(mat_det)
        return True, (rank, discriminant)

    return False, (rank, None)
    

if __name__ == "__main__":
    ############# An example of a conference Matrix #################
    ## We will create one that is 126 by 126 ##
    ## Notice that 126-1 = 5^3 so we can consider a degree 3 field extension of F_5
    ## here is some inputs to use
    ## where s=sqrt(n-1) under the field i think.
    ## maybe should check
    # p, d,  n,  s, k, sqrt(2d-1),        n-1
    # 89,63, 126,6, 18,11.180339887498949,5^3      Yes n-1=125=36=6^2=s^2  (mod 11)
    # 11,63, 126,13,9, 11.180339887498949,5^3      Yes n-1=125=4=2^2=13^2=s^2 (mod 11)
    # 11,63, 126,2, 42,11.180339887498949,5^3
    # 7, 27, 54, 5, 9, 7.280109889280518, 53^1     Yes n-1 = 53 = 4 = 25 = 5^2

    f125 = fieldmath.FieldExtension(fieldmath.Zp(5),[1,1,0,1])
    #f25 = fieldmath.FieldExtension(fieldmath.Zp(5),[3,0,1])
    f11 = fieldmath.Zp(11)
    f89 = fieldmath.Zp(89)

    ### Lets start by looking F_11
    ### we will construct a 63 X 126 (1,3,2)-ETF
    d = 63
    n = 126  # 2d
    G = create_gram_of_d_2d_etf_from_conference_mat(create_conference_matrix(f125, f11), f11, 2, True)
    #print(is_etf(None, G, True))
    # I still havent verified this is a frame. But I can by looking at if the discrimenant, which really comes down to if the gram matrix of the IP on the image is invertible
    # But it is
    #is_frame(gram_mat=G, with_discr=False)

    ## Now we can look for a regular s-simplex
    s = 13  # or also 2? 

    # We can test some possible frame vectors
    its = 0
    print((math.comb(G.rows, s+1)//1000))

    for sub_frame_G, sub_frame_cols in fieldmath.iter_sub_mats(G, s+1):
        if its % 1000 == 0:
            print(".", end='', flush=True)
        its+=1

        if sub_frame_G.rank() != s:
            continue

        is_tgt, c = is_tight(None, sub_frame_G)
        if is_tgt:
            print(f"Yo! {sub_frame_cols} is a {c}-tight frame, therefore a regular s-simplex!")
            break
    print("Sad stuff, no regular simplex here")

