import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math
import itertools
import numpy as np


bf = fieldmath.Zp(5)
f3 = fieldmath.Zp(3)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)

def find_incoherence(Phi=None, gram_mat=None):

    if Phi is not None:
        d = Phi.rows
        n = Phi.columns
        f = Phi.f
        gram_mat = Phi.adjoint()*Phi
    else:
        d = d # im trusting you here
        n = gram_mat.columns
        f = gram_mat.f

    (isit, (a,b,c)) = frame_thy.is_etf(Phi=Phi, gram_mat=gram_mat)
    assert isit, "uhh really? you think you could trick me like that!"

    sqrs = [(x,f.print_elm(f.multiply(x, x))) for x in f.iter_elems()]
    def square_root(x):
        return next((srx[0] for srx in sqrs if srx[1]==f.print_elm(x)), None)
    def is_nonzero_sqr(x):
        sqroot = square_root(x)
        if sqroot is None or f.equals(sqroot, f.zero()):
            return False
        else:
            return True

    PosTripleCounter = 0
    mat_list = []
    for ii in range(0, n-2):
        for jj in range (ii+1, n-1):
            for kk in range(jj+1, n):
                if is_nonzero_sqr(f.multiply(f.multiply(gram_mat.get(ii,jj),gram_mat.get(jj,kk)), gram_mat.get(kk,ii))):
                    PosTripleCounter += 1
                    row_lst = [0]*n
                    row_lst[ii] = 1
                    row_lst[jj] = 1
                    row_lst[kk] = 1
                    mat_list.append(row_lst)
        
    # Incrementally finding j-tuples for j>3
    jTuple = np.array(mat_list)
    for jTupleSize in range(3, d):
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

    return jTuple


#print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

famed_4_10_etf = fieldmath.create_matrix(
    [
        [0, 0, 0,  0,  1,  1,  1,  1,  1,  1],
        [0, 0, 1,  1,  0,  0,  1,  1,  2,  2],
        [1, 1, 0,  0,  0,  0,  1,  2,  1,  2],
        [1, 2, 1,  2,  1,  2,  0,  0,  0,  0]
    ], f3, column_space_gram=[1,1,1,2])

merc_frame = fieldmath.create_matrix(
    [
        [0, 6, 5],
        [2, 10, 10]
    ], fieldmath.Zp(11)
)

merc_frame_13 = fieldmath.create_matrix(
    [
        [0, -4, 4],
        [2, 12, 12]
    ], fieldmath.Zp(13)
)

print(frame_thy.is_frame(merc_frame, with_discr=True))
print(frame_thy.is_etf(merc_frame))

print(frame_thy.contains_simplex(2, merc_frame))

print(find_incoherence(merc_frame))

#print(frame_thy.contains_simplex(17-3, idk2))