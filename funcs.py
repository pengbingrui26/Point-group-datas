import numpy as np

# =============================================================
sigma0 = np.array( [ [ 1,0 ], [0,1] ] ).astype( complex )
sigma1 = np.array( [ [ 0,1 ], [1,0] ] ).astype( complex )
sigma2 = np.array( [ [ 0,-complex(0,1) ], [complex(0,1),0] ] ).astype( complex )
sigma3 = np.array( [ [ 1,0 ], [0,-1] ] ).astype( complex )
sigma = [ sigma0, sigma1, sigma2, sigma3 ]

# ==============================================================

def exp_sl( pm, nfold ) :
	assert pm in [1,-1]
	return np.exp( pm*complex(0,1)*np.pi/nfold )

# ==============================================================

def exp_su2( pm, n, nfold, m, acc = 4 ) : # exp( pm i sigma*n * pi/nfold )
	assert pm in [1,-1], pm
	n = n/np.linalg.norm( np.array(n) )
	phi = np.pi/nfold*m
	ndotsigma = n[0] * sigma1 + n[1] * sigma2 + n[2] * sigma3
	mat =  np.cos(phi) * sigma0 + pm * complex(0,1) * np.sin(phi) * ndotsigma
	#return np.round( mat, acc )
	return mat

# =======================================================================

def ge_su2( vec ):
	vec = vec/np.linalg.norm( np.array(vec) )
	return vec[0] * sigma1 + vec[1] * sigma2 + vec[2] * sigma3

# =======================================================================
'''
def dsum( A, B ) : # direct sum 
	nA, nB = len(A), len(B) 
	C = np.zeros( (nA+nB,nA+nB) ).astype('complex')
	for i in range(nA) :
		C[i][:nA] = A[i]
	for j in range(nA,nA+nB) :
		print( C[j][nA:] )
		print( B[j-nA] )  
		C[j][nA:] = B[j-nA]
	return C
'''

def dsum( mat_list ) : # direct sum 
	assert len(mat_list) >=2, len(mat_list) 
	n = len(mat_list)
	if n == 2 :
		A, B = mat_list[0], mat_list[1]
		nA, nB = len(A), len(B) 
		C = np.zeros( (nA+nB,nA+nB) ).astype('complex')
		for i in range(nA) :
			C[i][:nA] = A[i]
		for j in range(nA,nA+nB) :
			C[j][nA:] = B[j-nA]
		return np.array( C )
	else :
		A, B = mat_list[0], dsum( mat_list[1:] )
		nA, nB = len(A), len(B) 
		C = np.zeros( (nA+nB,nA+nB) ).astype('complex')
		for i in range(nA) :
			C[i][:nA] = A[i]
		for j in range(nA,nA+nB) :
			C[j][nA:] = B[j-nA]
		return np.array( C )

#print( dsum( [ sigma0, sigma0, sigma0 ] ) )
#exit()

# ===================================================================

def dpro( mat_list ) :
	assert len(mat_list) >= 2, len(mat_list)
	n = len(mat_list)
	if n == 2 :
		return np.kron( mat_list[0], mat_list[1] )
	else :
		return np.kron( mat_list[0], dpro(mat_list[1:]) )

#print(dpro( [ sigma0, sigma0 ] ))
#exit()

# ===================================================================
'''
M = np.kron( sigma0, complex(0,1)*sigma3 )
T = np.kron( complex(0,1)*sigma2, sigma0 )
#print( T )
u, v = np.linalg.eig( sigma2 )
#print( np.dot( v.T, np.dot( sigma2, np.linalg.inv( v.T ) ) )  )
V = np.kron( v, sigma0 )
TV = np.dot( V.T, np.dot( T, np.linalg.inv( V.T ) ) ) 
TV = np.round( TV, 6 )
print(TV)
MV = np.dot( V.T, np.dot( M, np.linalg.inv( V.T ) ) ) 
MV = np.round( MV, 6 )
print(MV)
'''

#print( exp_su2( [1,1,1], 3 ) ) 

'''
ma = np.array( [ [ np.exp( complex(0,1)*5/12*np.pi ), \
				np.exp( -complex(0,1)/12*np.pi ) ], \
				 [ np.exp( complex(0,1)*5/12*np.pi ), \
				np.exp( complex(0,1)*11/12*np.pi ) ] ] )
ma = ma*np.sqrt(2)/2
print( ma )
print( np.round( np.linalg.eig( ma )[0], 4 ) )
'''


'''
c3 = exp_su2( 1, [1,1,1], 3, 1, 1 )
#c3 = -sigma0
C3 = np.kron( [ [1,0,0,0], \
				[0,0,0,1],\
				[0,1,0,0],\
				[0,0,1,0] ], c3 )
#C3 = np.kron( sigma0, C3 )

c2 = exp_su2( 1, [0,0,1], 3, 1, 1 )
C2 = np.kron( [ [0,1,0,0],
				[1,0,0,0],
				[0,0,0,1],
				[0,0,1,0] ], c2 )
#C2 = np.kron( sigma0, C2 )


T = np.kron( sigma0, np.kron( sigma0, complex(0,1)*sigma2 )  )
#assert len(T) == len(C2)
#assert len(T) == len(C3)
#comm_C2_T = np.dot( C2, T ) - np.dot( T, C2.conjugate() ) 
#comm_C3_T = np.dot( C3, T ) - np.dot( T, C3.conjugate() ) 

#print( np.round( comm_C2_T, 4 ) )
#print( np.round( comm_C3_T, 4 ) )
#print( np.linalg.norm( comm_C2_T ) )
#print( np.linalg.norm( comm_C3_T ) )

for i in range(4) :
	for j in range(4):
		for k in range(4):
			S = np.kron( sigma[i], np.kron( sigma[j], sigma[k] ) )
			comm_C3 = np.dot( S, C3 ) - np.dot( C3, S )
			comm_C2 = np.dot( S, C2 ) - np.dot( C2, S )
			#comm_T = np.dot( tmp, T ) - np.dot( T, tmp.conjugate() )
			#if abs(np.linalg.norm(comm_T))< 1e-4 :
			#if abs(np.linalg.norm(comm_C3))< 1e-4 :
			if abs(np.linalg.norm(comm_C2))< 1e-4 :
				print( i,j,k )
'''



'''
C3 = np.kron( sigma0, C3 )
C2 = np.kron( sigma0, C2 )  		
'''

# select S to satisfy [ C3, S ] = [ C2, S ] = 0

'''
for i in range(3) :
	for j in range(3):
		for k in range(3):
			for l in range(3) :
				tmp = np.kron( sigma[i], np.kron( sigma[j], sigma[k] ) )
				tmp = np.kron( tmp, sigma[l] )
				comm_C2 = np.dot( tmp, C2 ) - np.dot( C2, tmp )
				comm_C3 = np.dot( tmp, C3 ) - np.dot( C3, tmp )
				if abs(np.linalg.norm(comm_C2))< 1e-4 :
					if abs(np.linalg.norm(comm_C3))< 1e-4 :
						print( i,j,k,l )
'''

'''
T = np.kron( sigma0, \
		np.kron( sigma0, \
			np.kron( sigma0, complex(0,1)*sigma2 )  ) ) 

#S = np.kron( sigma1, \
#		np.kron( sigma0, \
#			np.kron( sigma0, sigma0 )  ) ) 
S = np.kron( sigma2, \
		np.kron( sigma0, \
			np.kron( sigma0, sigma0 )  ) ) 
'''

'''
# select mass term 
for i in range(3) :
	for j in range(3):
		for k in range(3):
			for l in range(3) :
				tmp = np.kron( sigma[i], np.kron( sigma[j], sigma[k] ) )
				tmp = np.kron( tmp, sigma[l] )
				comm_C2 = np.dot( tmp, C2 ) - np.dot( C2, tmp )
				comm_C3 = np.dot( tmp, C3 ) - np.dot( C3, tmp )
				comm_S = np.dot( tmp, S ) + np.dot( S, tmp )
				comm_T = np.dot( tmp, T ) - np.dot( T, tmp.conjugate() )
				if abs(np.linalg.norm(comm_C2))< 1e-4 :
					if abs(np.linalg.norm(comm_C3))< 1e-4 :
						if abs(np.linalg.norm(comm_S))< 1e-4 :
							if abs(np.linalg.norm(comm_T))< 1e-4 :
								print( i,j,k,l )
'''



'''
#print( np.round(v2, 6 ) )
C2U = np.dot( u2.T, np.dot( C2, np.linalg.inv(u2.T) ) )
#print( np.round( C2U, 2 ) )
C3U = np.dot( u2.T, np.dot( C3, np.linalg.inv(u2.T) ) )
#print( np.round( C3 ) )  
#print( np.round( C3U, 6 ) )

#C3 = np.array( [ [0,0,1],[1,0,0],[0,1,0] ] )
#C2 = np.array( [ [0,1,0],[1,0,0],[0,0,1] ] )
#print( np.dot( C3, np.dot( C2, np.linalg.inv(C3) ) ) )
'''

# =======================================

'''
T = dpro( [ sigma0, sigma0, complex(0,1)*sigma2 ] )
#S4 = dsum( [ exp_su2( 1, [0,0,1], 4, 1, 1 ), \
#			 exp_su2( -1, [0,0,1], 4, 1, 1 ), \
#			 exp_su2( 1, [0,0,1], 4, 3, 1 ), \
#			 exp_su2( -1, [0,0,1], 4, 3, 1 )  ] )

S4 = np.kron( np.array( [ [0,1,0,0], \
						  [0,0,1,0], \
						  [0,0,0,1], \
						  [1,0,0,0] ] ), exp_su2( 1, [0,0,1], 4, 1, 1 ) )

print( np.round( np.linalg.eig(S4)[0] , 6 ) )

S4 = np.array(S4) 
#comm = np.dot( T, S4.conjugate() ) - np.dot( S4, T )
#print( comm )

for i in range(4) :
	for j in range(4) :
		for k in range(4) :
			tmp = np.kron( sigma[i], np.kron( sigma[j], sigma[k] ) ) 
			comm_S4 = np.dot( tmp, S4 ) + np.dot( S4, tmp )
			if abs( np.linalg.norm(comm_S4) ) <= 1e-4 :
				print( 'S:', i, j ,k )
				for a in range(4) :
					for b in range(4) :
						for c in range(4) :
							mass = np.kron( sigma[a], \
								np.kron( sigma[b], sigma[c] ) ) 
							com_S4 = np.dot( mass, S4 ) - np.dot( S4, tmp )
							com_T = np.dot( mass, T ) - np.dot( T, mass.conjugate() )
							com_tmp = np.dot( mass, tmp ) + np.dot( tmp, mass )
							if all( [ abs( np.linalg.norm( xx ) ) < 1e-4 for xx in \
										[ com_T, com_S4, com_tmp  ] ] ) :
								print( 'mass:', a, b, c )
				print( '\n' )					
'''	


'''
C3 = np.kron( [ [1,0,0,0], \
				[0,0,1,0], \
				[0,0,0,1], \
				[0,1,0,0] ], exp_su2( 1, [1,1,1], 3, 1, 1 ) )

#C3 = np.kron( [ [1,0,0,0], \
#				[0,0,1,0], \
#				[0,0,0,1], \
#				[0,1,0,0] ], -sigma0 )

S = dpro( [ sigma0, sigma0, exp_su2( 1, [1,1,1], 3, 1, 1 ) ] )

C2 = np.kron( [ [0,1,0,0], \
				[-1,0,0,0], \
				[0,0,0,1], \
				[0,0,-1,0] ], sigma0 )

print( C2 - np.kron( sigma0, np.kron( complex(0,1)*sigma2, sigma0 ) ) )
#print( np.linalg.eig(C2)[0] )

V, U = np.linalg.eig( C3 )

C3U = np.dot( np.linalg.inv(U), np.dot( C3, U ) )
#print(np.round(C3U,1))
C2U = np.dot( np.linalg.inv(U), np.dot( C2, U ) )
#print( np.round(C2U,2) )
SU = np.dot( np.linalg.inv(U), np.dot( S, U ) )

comm_S_C3 = np.dot( S, C3 ) - np.dot( C3, S )
#print( comm_S_C3 )
#print( np.round(SU,1) )

#exit()

for i in range(4):
	for j in range(4):
		for k in range(4):
			S = dpro( [ sigma[i], sigma[j], sigma[k] ] )
			S = np.array(S)
			assert S.shape == C3.shape
			comm_C3 = np.dot( S, C3 ) - np.dot( C3, S )
			#if abs(np.linalg.norm(comm_C3)) < 1e-4 :
			#	print( i,j,k ) 
			comm_C2 = np.dot( S, C2 ) - np.dot( C2, S )
			if abs(np.linalg.norm(comm_C2)) < 1e-4 :
				print( i,j,k ) 
'''

# ===================================================================

'''
C3_111 = dpro( [ np.array( [ [1,0,0,0], [0,0,1,0],[0,0,0,1],[0,1,0,0] ] ), \
				exp_su2( 1, [1,1,1], 3, 1) ] )
#print( C3_111 )

s1 = ge_su2( [1,1,1] )
s2 = ge_su2( [-1,-1,1] )
s3 = ge_su2( [1,-1,-1] )
s4 = ge_su2( [-1,1,-1] )
S = dsum( [ s1, s2, s3, s4 ] ) 

T = dpro( [ sigma0, sigma0, complex(0,1)*sigma2 ] )

#v, u = np.linalg.eig(S)
#print( v )
print( np.round( np.dot( C3_111, S ) - np.dot( S, C3_111 ), 2 ) )

V, U = np.linalg.eig( C3_111 )
C3_111U = np.dot( np.linalg.inv(U), np.dot( C3_111, U ) )
print( np.round( C3_111U, 1 ) )
SU = np.dot( np.linalg.inv(U), np.dot( S, U ) )
VV, UU = np.linalg.eig( SU )
SUU = np.dot( UU.T, np.dot( SU, np.linalg.inv(UU.T) ) )
print( np.round( SUU, 1 ) )
C3_111UU = np.dot( UU.T, np.dot( C3_111U, np.linalg.inv(UU.T) ) )
print( np.round( C3_111U, 1 ) )
'''
