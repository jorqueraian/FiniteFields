import fieldmath
import itertools
import math
import numpy as np


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

    if field.equals(sqrt_2dminus1, 0):
        return conference_mat + fieldmath.identity_n(conference_mat.rows, field, sqrt_2dminus1)
    else:
        return conference_mat*field.reciprocal(sqrt_2dminus1)+ fieldmath.identity_n(conference_mat.rows, field)


def is_equiangular(Phi=None, gram_mat=None, why_not=False):
    """Assumes equinorm, and then computes (a,b) parameters"""
    if Phi is None and gram_mat is None:
        assert False, "you need to specify an input"
    
    if gram_mat is None:
        gram_mat = (Phi.adjoint()*Phi)
    gram_mat_modulus_sqrd = gram_mat.modulus_squared_of_entries()
    a = gram_mat.get(0, 0) 
    a_sqrd = gram_mat_modulus_sqrd.get(0, 0)
    b = gram_mat_modulus_sqrd.get(0, 1)

    # This is lazy but im ok with that.
    # I have changed a, did this break things?
    test_mat = fieldmath.identity_n(gram_mat.rows, gram_mat.f, gram_mat.f.subtract(a_sqrd,b)) + fieldmath.Matrix(gram_mat.rows, gram_mat.rows, gram_mat.f, b)
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
        gram_mat = (Phi.adjoint()*Phi)

    gram_sqrd = gram_mat*gram_mat

    if gram_mat.f.equals(gram_sqrd.get(0,0), gram_mat.f.zero()) and gram_sqrd.any():
        if why_not:
            print("Yikes")
        return False, None
    elif not gram_sqrd.any():
        return True, gram_sqrd.get(0,0)
    
    if gram_mat.f.equals(gram_mat.get(0,0), gram_mat.f.zero()):
        return False, None
    c = gram_mat.f.divide(gram_sqrd.get(0,0),gram_mat.get(0,0))

    if (gram_sqrd-(gram_mat*c)).any():
        if why_not:
            print("Its just not in the cards, sorry.")
        return False, None
    else:
        return True, c
    
def spark(Phi):
    # This should only every be used for small examples. This is a NP-hard problem
    # so brute force is basically all we got.
    d = Phi.rows
    n = Phi.columns
    n_list = [i for i in range(n)]
    spark = d+1
    done = False
    for num_to_test in range(d, 1, -1):
        if done == True:
            break
        done = True
        for subset in itertools.combinations(n_list, num_to_test):
            rnk = (Phi.get_sub_matrix_from_cols(subset)).rank()
            if rnk < len(subset):
                # This measn LD
                spark -= 1
                done = False
                break

    return spark


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
    
    if Phi is not None:
        fld = Phi.f
        rank = Phi.rank()
        if gram_mat is None:
            gram_mat = (Phi.adjoint()*Phi)
        frame_bool = (Phi.rank() == (Phi.adjoint()*Phi).rank())
    elif gram_mat is not None:
        fld = gram_mat.f

        if (gram_mat-gram_mat.adjoint()).any():
            return False, (None, None)
        frame_bool = True
        rank = gram_mat.rank()

    if with_discr is False or not frame_bool:
        return frame_bool, (rank if frame_bool else None, None)
    else:
        sqrs = list(set([fld.print_elm(fld.multiply(x, x)) for x in fld.iter_elems()]))
        
        def is_sqr(x):
            if fld.equals(x, fld.zero()):
                return fld.zero()
            if (fld.print_elm(x) in sqrs):
                return True
            else:
                return False

        # only need to do this if the field automorphism is trivial
        discriminant = None
        total_columns_set = {k for k in range(gram_mat.columns)}

        # Find basic submatrix.
        print("THIS IS SLOW, probably")
        for columns_lst in itertools.combinations(total_columns_set, rank):
            selected_columns = gram_mat.get_sub_matrix_from_lists(columns_lst, columns_lst)
            if selected_columns.rank() != rank:
                continue

            sub_mat = gram_mat.get_sub_matrix_from_lists(columns_lst, columns_lst)
            mat_det = sub_mat.det()
            discriminant = is_sqr(mat_det)
            return True, (rank, discriminant)
    

def contains_simplex(s, Phi=None, gram_mat=None):
    """function Binder = BinderFinder(Phi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Binder Finder
    %%%%%% Code to find the embedded simplices (ETFs for their spans with one
    %%%%%% more vector than their rank) in a given ETF.
    %%%%%%
    %%%%%% Input: Phi, the synthesis matrix of an equiangular tight frame.
    %%%%%%
    %%%%%% Output: Binder, the incidence matrix (rows = simplices, cols = frame
    %%%%%% vectors) of the simplices embedded in the given ETF
    %%%%%%
    %%%%%% Citation: "Equiangular tight frames that contain regular simplices,"
    %%%%%% Matthew Fickus, John Jasper, Emily J. King, Dustin G. Mixon, 2017
    %%%%%%
    %%%%%% Code created: July 2016
    %%%%%% Last updated: November 22, 2017
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"""


    if Phi is not None:
        d = Phi.rows
        n = Phi.columns
        f = Phi.f
        gram_mat = Phi.adjoint()*Phi
    else:
        d = d # im trusting you here
        n = gram_mat.columns
        f = gram_mat.f
    
    k = s+1


    (isit, (a,b,c)) = is_etf(Phi=Phi, gram_mat=gram_mat)
    assert isit, "uhh really? you think you could trick me like that!"

    assert f.equals(f.multiply(a,a),f.multiply(f.multiply(s,s),b)), "Not even compatible, get out of here with those numbers"

    sqrs = [(x,f.print_elm(f.multiply(x, x))) for x in f.iter_elems()]
    def square_root(x):
        return next((srx[0] for srx in sqrs if srx[1]==f.print_elm(x)), None)
    
    if not f.equals(a, f.zero()):
        assert square_root(f.divide(a,s)) is not None, "OK really?"

        neg_a3bys3 = f.negate(fieldmath.pow_over_field(f.divide(a,s), 3, f, False))

        TripleCounter = 0
        mat_list = []
        for ii in range(0, n-2):
            for jj in range (ii+1, n-1):
                for kk in range(jj+1, n):
                    if f.equals(f.multiply(f.multiply(gram_mat.get(ii,jj),gram_mat.get(jj,kk)), gram_mat.get(kk,ii)), neg_a3bys3):
                        TripleCounter += 1
                        row_lst = [0]*n
                        row_lst[ii] = 1
                        row_lst[jj] = 1
                        row_lst[kk] = 1
                        mat_list.append(row_lst)
        
        # Incrementally finding j-tuples for j>3
        jTuple = np.array(mat_list)
        for jTupleSize in range(3, k):
            for ii in range(0, jTuple.shape[0]):  # should get number of rows or 
                Indices = np.where(jTuple[ii, :] == 1)[0]
                Indicator = np.matmul(np.transpose(np.sum((np.matmul(jTuple[:,Indices],(np.ones(jTupleSize)-np.eye(jTupleSize))) == jTupleSize-1).astype(int),axis=1)),jTuple)
                Indicator[Indices] = 0
                NewIndices = np.where(Indicator == jTupleSize)[0]
                A = np.kron(np.ones((NewIndices.shape[0],1),dtype=int),jTuple[ii,:])
                A[:,NewIndices] = np.eye(NewIndices.shape[0])
                if ii == 0:
                    NewjTuple = A
                else:
                    NewjTuple = np.vstack([NewjTuple, A])
                    #NewjTuple[(NewjTuple.shape[0]-1):(NewjTuple.shape[0]+A.shape[0]-1),:] = A
                jTuple[ii,:] = np.zeros((1,jTuple.shape[1]))

            jTuple = NewjTuple

        Binder = jTuple
    else:
        Binder = []
        for maybe_simplex_ind in itertools.combinations(range(Phi.columns), s+1):
            maybe_simplex = Phi.get_sub_matrix_from_cols(maybe_simplex_ind)
            simpl_gram = (maybe_simplex.adjoint()*maybe_simplex)
        
            if maybe_simplex.rank() != s or simpl_gram.rank() != maybe_simplex.rank():
                continue
                        
            # Now we know its a frame
        
            (is_c_simp, c_simp) = is_tight(gram_mat=simpl_gram)
                            
            if not is_c_simp:
                continue
            
            Binder.append(maybe_simplex_ind)
    return Binder#, all_simps


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

    f49 = fieldmath.FieldExtension(fieldmath.Zp(7),[1,0,3])

    ### Lets start by looking F_11
    ### we will construct a 63 X 126 (1,3,2)-ETF
    d = 25
    n = 50  # 50 = 7^2+1
    G = create_gram_of_d_2d_etf_from_conference_mat(create_conference_matrix(f49, fieldmath.Zp(7)), fieldmath.Zp(7), 0, True)
    #print(is_etf(None, G, True))
    # I still havent verified this is a frame. But I can by looking at if the discrimenant, which really comes down to if the gram matrix of the IP on the image is invertible
    # But it is
    #is_frame(gram_mat=G, with_discr=False)

    ## Now we can look for a regular s-simplex
    s = 21  # or also 21? 

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

