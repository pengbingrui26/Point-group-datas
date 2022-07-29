import numpy as np
import sympy as sp
from sympy.core.numbers import I
from functools import reduce
from scipy.linalg import null_space
import os
import itertools

from funcs import sigma0, sigma1, sigma2, sigma3, sigma, exp_sl, exp_su2, ge_su2, dsum, dpro

def float_to_sympy(a):
    # convert a float number to a sympy number
    def is_int(num):
        # return True if the input number is close enough to an integer
        return True if abs(num - round(num)) < 1e-8 else False

    if abs(a) < 1e-8:
        return 0
    sign = 1 if abs(abs(a) - a) < 1e-8 else -1 # use abs(diff)<1e-8: when a is a sp number, abs(a)==a is always False 
    a = abs(a)
    sqrt_dict = {1:1, 2:1.4142135624, 3:1.7320508075, 4:2, 5:2.23606797749, 6:2.449489742783178}
    denominator_num = [2,3,4,5,6]
    for num, sqrt_num in sqrt_dict.items():
        if is_int(sqrt_num / a):
            time = int(round(sqrt_num / a))
            return sp.sqrt(num) / time * sign
        elif is_int(a / sqrt_num):
            time = int(round(a / sqrt_num))
            return sp.sqrt(num) * time * sign
        
    for num, sqrt_num in sqrt_dict.items():
        for denom in denominator_num:
            if is_int(sqrt_num / a * denom):
                time = int(round(sqrt_num / a * denom))
                return sp.sqrt(num) / time * denom * sign
            elif is_int(a / sqrt_num * denom):
                time = int(round(a / sqrt_num * denom))
                return sp.sqrt(num) * time / denom *sign

    for num in frac_num:
        if is_int(a * num):
            time = int(round(a * num))
            return sp.Rational(time / num) * sign

    if abs(a - 0.6830127019) < 1e-8:
        return (sp.sqrt(3) + 1)/4 * sign
    elif abs(a - 0.1830127019) < 1e-8:
        return (sp.sqrt(3) - 1)/4 * sign
    else:
        raise ValueError('Fail to identify float,a='+str(a))

def numpy_to_sympy(A):
    # convert a numpy square matrix to sympy matrix
    shape = np.shape(A)
    return np.array([[ float_to_sympy(np.real(A[i,j])) + I * float_to_sympy(np.imag(A[i,j])) for j in range(shape[1]) ] for i in range(shape[0]) ])

def round_mat(A):
    # round the small decimal part of matrix A
    def round_num(a):
        # round the small decimal part of float a 
        real, imag = np.real(a), np.imag(a)
        real = np.round(real) if abs(real - np.round(real)) < 1e-8 else real
        imag = np.round(imag) if abs(imag - np.round(imag)) < 1e-8 else imag
        return real + imag * 1j

    shape = np.shape(A)
    return np.array([[ round_num(A[i,j]) for j in range(shape[1]) ] for i in range(shape[0]) ])

def mat_eq(M1, M2):
    return True if np.max(np.abs(M1 - M2)) < 1e-8 else False

def phase_correction(C):
    # try to remove the phase ambiguity of C and turn C into an integer matrix if possible
    for i in range(np.shape(C)[0]):
        for j in range(np.shape(C)[1]):
            if abs(C[i,j]) < 1e-8: #c=0
                C[i,j] = 0
                continue
            if abs(abs(C[i,j]) - 1) < 1e-8: #c=exp(i*theta)
                C = C/C[i,j]
    return round_mat(C)
    

def find_similarity_transformation(A,B):
    # find a unitary matrix C, s.t A = CBC^-1 ==> AC=CB
    # Method: tranform into a linear eq of C, and then solve C
    """
    Example:
    >>> A = [np.diag((1,-1)), np.array([[0,1],[1,0]])]
    >>> B = [np.diag((1,-1)),np.array([[0,1j],[-1j,0]])]
    >>> find_similarity_transformation(A,B)
    [[1.+0.j 0.+0.j],[0.+0.j 0.+1.j]]
    """
    def transform_AC_eq_CB(A,B):
        # transform AC-CB=0 into a linear eq of C, i.e, mat*C'=0, where C' is the flattened column vector of C.
        dim = np.shape(A[0])[0]
        num_mat = len(A)
        mat = np.zeros((num_mat * dim**2, dim**2), dtype=complex)
        for imat in range(num_mat):
            a, b = A[imat], B[imat]
            for i1 in range(dim):
                for i2 in range(dim):
                    for j1 in range(dim):
                        for j2 in range(dim):
                            if j2 == i2:
                                mat[imat*dim**2 + i1*dim+i2, j1*dim+j2] += a[i1, j1]
                            if j1 == i1:
                                mat[imat*dim**2 + i1*dim+i2, j1*dim+j2] -= b[j2,i2]
        return mat

    def transform_colC_to_matC(colC, dim):
        # transfrom N^2 dim column vec C to N*N unitary matrix C
        # shape(C) =(dim^2,m), where m is the number of indepent null vec
        assert np.shape(colC)[0] == dim**2 and np.shape(colC)[1] > 0, ('Wrong colC!',colC)
        trial_C_list = [ colC[:,ii].reshape((dim,dim)) for ii in range(np.shape(colC)[1])]
        trial_C_list.append(reduce(lambda x, y : x + y, trial_C_list))
        trial_C_list.append((colC[:,0] + colC[:,-1]).reshape((dim,dim)))

        for tmp in trial_C_list:
            CCdagger = np.dot(tmp, np.conjugate(tmp).T)
            if CCdagger[0,0].real != 0 and np.linalg.norm(CCdagger/CCdagger[0,0].real - np.eye(dim)) < 1e-6:
                return tmp / np.sqrt(CCdagger[0,0].real)
        else:
            raise ValueError('Fail to find unitary C! null vec='+str(colC))

    dim = np.shape(A[0])[0]
    trans_mat = transform_AC_eq_CB(A, B)
    tmpC = null_space(trans_mat)
    C = transform_colC_to_matC(tmpC, dim)
    C = phase_correction(phase_correction(C))
    assert sum([np.linalg.norm(np.dot(A[i],C)-np.dot(C,B[i]))>1e-6 for i in range(len(A))])<1e-6,'Wrong C! C='+str(C)
    return C

def rep_similarity_transformation(A, B, block_dim = None, nir_list = None, irlabel_list = None):
    """
    for two equivalent representations A, B, find unitary C s.t C^-1AC = B
    If B is irreducible, block_dim=None, while if B is reducible, block_dim is a list to denote the dim of each block irrep.
    Algorithm:
    I. If B is irreducible: 
         M := sum(A_i \otimes B_i^*)/g, where g=len(B)
         Calculate the eigenvector v with eigenvalue 1 of M (assert there exists one such v, and only one)
         C = v.reshape(n,n)/sqrt(g), where (n,n) is the shape of B
    II. If B is reducible, each B_i=(B_i^a1, B_i^a2, ...) must be block-diagonal, with the dim of each block given in block_dim=[d1,d2..]
         M^a := sum(A_i \otimes B_i^a*)/g, calculate eigenvector v^a with eigenvalue 1 of M^a
         C^a = v.reshape(n,d^a), 
         C=(C^a1, C^a2, ...)

    Limitations: This algorithm takes only rep matrix of unitary symmetries.
    For anti-unitary symmetries, the transformation rule is C^-1 A C^*=B, which have not been considered.
    """
    def eigvec_of_eigval_1(M):
        eig_vals, eig_vecs = np.linalg.eig(M)
        locate_1 = [ np.abs(d - 1) < 1e-8 for d in eig_vals ]
        assert sum(locate_1) > 0, ('M does not have eigenvalue 1', locate_1,eig_vals)
        if sum(locate_1) == 1:
            v = eig_vecs[:, locate_1.index(1)]
            return 1, [v]
        else:
            num_1 = sum(locate_1)
            v_list = [ eig_vecs[:,icol] for icol, col in enumerate(locate_1) if col == True ]
            return num_1, v_list

    assert len(A) == len(B) and np.shape(A[0]) == np.shape(B[0])
    nop, dim = len(A), np.shape(A[0])[0]
    if not block_dim:
        M = sum([ np.kron(Ai, np.conjugate(Bi)) for Ai, Bi in zip(A, B) ]) / nop
        num_1, v = eigvec_of_eigval_1(M)
        assert num_1 == 1
        C = v[0].reshape((dim, dim)) * np.sqrt(dim)
    else:
        C = np.zeros((dim, dim),dtype=complex)
        for ith, IRdim, IRlabel in zip(range(len(block_dim)), block_dim, irlabel_list):
            start = sum([ d * n for d, n in zip(block_dim[:ith], nir_list[:ith]) ])
            if '_' not in IRlabel:
                end = start + IRdim
                Ma = sum([ np.kron(Ai, np.conjugate(Bi[start:end, start:end])) for Ai, Bi in zip(A, B) ]) / nop
                num_1, va = eigvec_of_eigval_1(Ma)
                assert num_1 == nir_list[ith], (num_1, nir_list, ith, IRlabel)
                for nth, v in enumerate(va):
                    Ca = v.reshape((dim, IRdim))
                    C[:, start + nth * IRdim : end + nth * IRdim] = Ca * np.sqrt(IRdim)

            else: # coir of two irreps
                assert IRdim % 2 == 0, ('In this case, IRdim should be even', IRdim)
                subdim = int(IRdim/2)
                if IRlabel.split('_')[0] != IRlabel.split('_')[1]: # two different ir
                    for icoir in range(2):
                        start = start + icoir * subdim
                        end = start + subdim
                        Ma = sum([ np.kron(Ai, np.conjugate(Bi[start:end, start:end])) for Ai, Bi in zip(A, B) ]) / nop
                        num_1, va = eigvec_of_eigval_1(Ma)
                        assert num_1 == nir_list[ith], (num_1, nir_list, ith)
                        for nth, v in enumerate(va):
                            Ca = v.reshape((dim, subdim))
                            C[:, start + nth * IRdim : end + nth * IRdim] = Ca * np.sqrt(subdim)
                else:
                    end = start + subdim
                    Ma = sum([ np.kron(Ai, np.conjugate(Bi[start:end, start:end])) for Ai, Bi in zip(A, B) ]) / nop
                    num_1, va = eigvec_of_eigval_1(Ma)
                    assert num_1 == 2 * nir_list[ith], (num_1, nir_list, ith)
                    for nth, v in enumerate(va):
                        Ca = v.reshape((dim, subdim))
                        C[:, start + nth * subdim : end + nth * subdim] = Ca * np.sqrt(subdim)

    C = phase_correction(C)
    assert mat_eq(C @ np.conjugate(C).T, np.eye(dim)), ('C not unitary!', C, C @ np.conjugate(C).T)
    assert all([ mat_eq(Ai @ C, C @ Bi) for Ai, Bi in zip(A, B) ]), 'C fails to tranform A to B!'        
    return C  

# ===============================================================   

#c3 = np.array( [ [1,0,0,0], [0,0,1,0],[0,0,0,1],[0,1,0,0] ] )
#s = np.eye( [] )

C3 = dpro( [ np.array( [ [1,0,0,0], [0,0,1,0],[0,0,0,1],[0,1,0,0] ] ), \
                 exp_su2( 1, [1,1,1], 3, 1) ] )
 
s1 = ge_su2( [1,1,1] )
s2 = ge_su2( [-1,-1,1] )
s3 = ge_su2( [1,-1,-1] )
s4 = ge_su2( [-1,1,-1] )
S = dsum( [ s1, s2, s3, s4 ] ) 

VC3, UC3 = np.linalg.eig(C3)
VS, US = np.linalg.eig(S)

#print(VC3[4])
#print(np.round(VC3,6))
#print(np.round(VS,6))

VS = np.round( VS, 6 )

pers = list( itertools.permutations(VS, len(VS)) )

pos = []

for xx in pers:
	if list(xx) not in pos:
		pos.append( list(xx) )		

C3_new = np.zeros( (8, 8) ).astype('complex')
for i in range(8):
	C3_new[i][i] = VC3[i]

#print( C3_new )

for k in range(len(pos)):		
	print( 'k:', k )
	S_new = np.zeros( (8, 8) ).astype('complex')	
	for ii in range(8):		
		S_new[ii][ii] = pos[k][ii]		
	print( S_new )
	commu = np.dot( C3_new, S_new ) - np.dot( S_new, C3_new )
	print( np.round( np.linalg.norm(commu), 6 ) )
	#A = [ C3, S ]
	#B = [ C3_new, S_new ]
	#B= A
#try:
#C = find_similarity_transformation(A,B)
	#C = rep_similarity_transformation(A,B)
#except:
#	continue
	#print( 'k:', k )
	#print(C)
	print('\n')

exit()

V, U = np.linalg.eig( S )
C3U = np.dot( np.linalg.inv(U), np.dot( C3, U ) )
print( np.round( C3U, 2 ) )
SU = np.dot( np.linalg.inv(U), np.dot( S, U ) )
print( np.round( SU, 2 ) )

#VV, UU = np.linalg.eig( SU )
#SUU = np.dot( UU.T, np.dot( SU, np.linalg.inv(UU.T) ) )
#print( np.round( SUU, 1 ) )
#C3_111UU = np.dot( UU.T, np.dot( C3_111U, np.linalg.inv(UU.T) ) )
#print( np.round( C3_111U, 1 ) )

