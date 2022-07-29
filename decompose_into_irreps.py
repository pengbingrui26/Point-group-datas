import numpy as np
from all_pg_irreps import all_pg_irreps
#from all_1d_irreps import all_1d_irreps
from irreps import getPG
from funcs import exp_su2

sigma0 = np.array( [ [1,0] ,[0,1] ] )
sigma1 = np.array( [ [0,1] ,[1,0] ] )
sigma2 = np.array( [ [0,-complex(0,1)] ,[complex(0,1),0] ] )
sigma3 = np.array( [ [1,0] ,[0,1] ] )
sigma = [ sigma0, sigma1, sigma2, sigma3 ]

# =======================================================================

c3 = exp_su2( 1, [1,1,1], 3, 1 )
#C3 = np.kron( np.array( [ [1,0,0,0],[0,0,1,0],[0,0,0,1],[0,1,0,0] ] ), c3 )
C3 = np.kron( np.array( [ [1,0,0,0],[0,0,1,0],[0,0,0,1],[0,1,0,0] ] ), \
				-sigma0 )

c2 = np.array( [ [ np.exp(complex(0,1)*np.pi/2), 0 ], \
					[ 0, np.exp(-complex(0,1)*np.pi/2) ] ] )
#C2 = np.kron( np.array( [ [0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0] ] ), c2 )
C2 = np.kron( np.array( [ [0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0] ] ), sigma0 )

C3inv = np.linalg.inv(C3)

Id = np.kron( sigma0, np.kron( sigma0, sigma0 ) )


pg_T = all_pg_irreps[27]
#print( len((pg_T )['sym_ops']) ) 
#print( pg_T['dv_irrep_names'] )
nop = len(pg_T['sym_ops'])
nirrep = len(pg_T['dv_irrep_names'])
irreps = [ [ pg_T['dv_irreps'][i][k] for i in range(nop) ] for k in range(nirrep) ]
#print( len(irreps) )
#print( len(irreps[0]) ) 
PG = getPG(28)
#print( PG.cls ) 
cls_reps = [ cls[0] for cls in PG.cls ]
#print( cls_reps )
#for x in cls_reps :	
#	print( pg_T['sym_ops'][x][1] )
for i, irrep in enumerate(irreps) :
	print( pg_T['dv_irrep_names'][i] )
	irep_Id, irep_C2, irep_C3, irep_C3inv = [ irrep[j] for j in cls_reps ]
	#print( irep_C3 )  
	tot = 0
	tot += ( 1*np.trace(irep_Id.conjugate()) * np.trace(Id)   \
			+ 3*np.trace(irep_C2.conjugate()) * np.trace(C2)   \
			+ 4*np.trace(irep_C3.conjugate()) * np.trace(C3)  \
			+ 4*np.trace(irep_C3inv.conjugate()) * np.trace(C3inv)  )  
	tot/= 12
	real = tot.real
	imag = tot.imag
	assert abs( real-round(real) ) < 1e-4
	assert abs( imag ) < 1e-4
	tot = int(round(real))
	print( tot )

# ===========================================================================

'''
Mx = np.kron( sigma0, complex(0,1)*sigma1 )
My = np.kron( sigma1, complex(0,1)*sigma2 )
C2 = np.dot( Mx, My ) 
V,U = np.linalg.eig(C2)
MxU = np.dot( np.linalg.inv(U), np.dot(Mx,U) )
MyU = np.dot( np.linalg.inv(U), np.dot(My,U) )
C2U = np.dot( np.linalg.inv(U), np.dot(C2,U) )
print(np.round(MxU,2))
print(np.round(MyU,2))
print(np.round(C2U,2))
'''

# ===========================================================================

'''
def decompose_rep( irrep, rep ) :
	# n_irrep = 1/NG * ( \sum_G X_irrep^{*}(G) * X_rep(G) )
	assert len(irrep) == len(rep) 
	NG = len(irrep)
	tot = 0
	for i in range(NG) :
		tot += np.trace( irrep[i].conjugate() ) * np.trace( rep[i] )
	real = tot.real
	imag = tot.imag
	assert abs(imag) < 1e-4, tot
	assert abs( real - round(real) ) < 1e-4
	tot = int(round(real))
	tot /= NG
	return tot 
'''	 


