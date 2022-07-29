import numpy as np
from PG import symop, Pointgroup, getPG 
#from funcs import exp_su2
# =============================================================

sigma0 = np.array( [ [ 1,0 ], [0,1] ] ).astype( complex )
sigma1 = np.array( [ [ 0,1 ], [1,0] ] ).astype( complex )
sigma2 = np.array( [ [ 0,-complex(0,1) ], [complex(0,1),0] ] \
                         ).astype( complex )
sigma3 = np.array( [ [ 1,0 ], [0,-1] ] ).astype( complex )
sigma = [ sigma0, sigma1, sigma2, sigma3 ]

# ==============================================================

# ===============================================================================

def exp_su2( pm, n, nfold, m, acc = 4 ) : # exp( pm i sigma*n * pi/nfold )		
	assert pm in [1,-1], pm	
	n = n/np.linalg.norm( np.array(n) )		
	phi = np.pi/nfold*m		
	ndotsigma = n[0] * sigma1 + n[1] * sigma2 + n[2] * sigma3		
	mat =  np.cos(phi) * sigma0 + pm * complex(0,1) * np.sin(phi) * ndotsigma		
	#mat =  np.cos(phi) * sigma0 - pm * complex(0,1) * np.sin(phi) * ndotsigma		
	#return np.round( mat, acc )		
	return mat

# ===============================================================================

def ebr( PG, site_gen, site_irrep ) :
	# site_gen: generators for site_symmetry point group of the given Wykoff position
	# site_irrep: the irrep of the given Wykoff position
	assert type(PG) == Pointgroup
	#
	#for iop, op in enumerate(PG.op) :
	#	print( iop, op.name )
	#
	dim_irrep = len( site_irrep.values()[0] )   # dimension of the given irrep
	site_group = PG.get_subgroup( site_gen )
	coset_dict = PG.get_coset( site_group )
	cosets = [  ]
	for coset_rep in coset_dict :
		cosets.append( ( coset_rep, coset_dict[coset_rep] ) ) 
	#print( coset_dict )
	n_orbit = len(cosets)  # number of orbitals 
	#exit()

	for iop, op in enumerate(PG.op) :  # for each symop of PG, decide its representation matrix
		#if iop != 4 :  # C_{3,111,+} for T
		#if iop != 1 :  # C_{2z} for T 
		#if iop != 3 : # M_x for C2v
		#if iop != 2 : # C_2 for S4
		#	continue 
		print( 'op:', op.name )
		mat_orbit = np.zeros( (n_orbit,n_orbit) ).astype('complex')
		mat_ebr = np.zeros( (n_orbit*dim_irrep,n_orbit*dim_irrep) ).astype('complex')
		# for each orbit, decide which one it is transformed into by op 
		for i, coset in enumerate( cosets ) : 
			##coset_new = [ PG.mut[ ( xx, iop ) ] for xx in coset[1] ]
			coset_new = [ PG.mut[ ( iop, xx ) ] for xx in coset[1] ]
			coset_new.sort()
			i_new = 0
			find = False
			for j, coset_tmp in enumerate( cosets ) :
				if coset_tmp[1] == coset_new :
					i_new = j
					find = True
					break
			assert find
			g_alpha = coset[0]
			g_beta = cosets[i_new][0]
			print( 'g_alpha, g_beta:', g_alpha, g_beta )
			#print( 'g_beta^{-1} * h * g_beta:', \
			#		PG.op[PG.inv[g_beta]].name, PG.op[iop].name, PG.op[g_alpha].name )
			
			#g_su2 = np.dot( PG.op[PG.inv[g_beta]].su2, \
			#			np.dot( PG.op[iop].su2, PG.op[g_alpha].su2 ) )
			g_su2 = np.dot( np.linalg.inv(PG.op[g_beta].su2), \
						np.dot( PG.op[iop].su2, PG.op[g_alpha].su2 ) )
			gg_su2 = np.dot( PG.op[PG.inv[g_beta]].su2, \
						np.dot( PG.op[iop].su2, PG.op[g_alpha].su2 ) )
			ggg_su2 = PG.op[ PG.mut[ ( PG.inv[g_beta], PG.mut[ ( iop, g_alpha ) ] ) ] ].su2
			ggg = PG.mut[ ( PG.inv[g_beta], PG.mut[ ( iop, g_alpha ) ] ) ]
			#
			g_iop = ''
			find_g = False
			su2_pm = ''  
			for iii, oppp in enumerate(PG.op) :
				if abs(np.linalg.norm(oppp.su2-g_su2)) < 1e-4 :
					g_iop = iii
					find_g = True
					su2_pm = 1
					break
				elif abs(np.linalg.norm(oppp.su2+g_su2)) < 1e-4 :
					g_iop = iii
					find_g = True
					su2_pm = -1
					break
			assert find_g
			#
			print( 'g_su2:', np.round( g_su2 ) )
			print( 'gg_su2:', np.round( gg_su2 ) )
			print( 'ggg_su2:', np.round( ggg_su2, 2 ) )
			print( 'ggg:', ggg )
			# 
			assert g_iop in site_irrep.keys(), ( g_iop, site_irrep.keys() )
			irrep = site_irrep[g_iop] * su2_pm
			#
			mat_orbit[i][i_new] = 1
			for ii in range(dim_irrep) :
				for jj in range(dim_irrep) :
					mat_ebr[i*dim_irrep+ii][i_new*dim_irrep+jj] = irrep[ii][jj]
		print( 'mat_orbit:' )
		print(  mat_orbit )
		print( 'mat_ebr:' )
		#print( mat_ebr )
		print( np.round( mat_ebr, 6 ) )
		print( 'mat_ebr eigenvalues:' )
		print( np.linalg.eig( mat_ebr )[0] )
		print( '\n' )

# =====================================================================
# test EBR
# ====================================================================

# point group T 

'''
T = getPG(28)
site_gen = [4]
#ebr( T, site_gen )
for iop, op in enumerate(T.op):
	print( iop, op.name )

C3111p = exp_su2( 1, [1,1,1], 3, 1 ) 
C3111m = np.linalg.inv( C3111p )
Id = sigma0
#Idm = -sigma0

site_irrep = { 0:Id, 4:C3111p, 8:C3111m }   # { iop: irrep_mat, ....... } 
print( 'mtb:' )
print( T.mut )
ebr( T, site_gen, site_irrep )
'''

# point group C2v

'''
C2v = getPG(7)
site_gen = [3]
#ebr( T, site_gen )
for iop, op in enumerate(C2v.op):
	print( iop, op.name )

Mx = exp_su2( 1, [1,0,0], 2, 1 ) 
Id = sigma0

site_irrep = { 0:Id, 3:Mx }   # { iop: irrep_mat, ....... } 
ebr( C2v, site_gen, site_irrep )
'''

# point group S4

'''
S4 = getPG(10)
print( S4.inv )
print( S4.cls )
#site_gen = [1] # C2
site_gen = [0] # Id
for iop, op in enumerate( S4.op ) :
	print( iop, op.name )

C2 = exp_su2( 1, [0,0,1], 2, 1 )
Id = sigma0

#site_irrep = { 0: Id, 1: C2 }
#ebr( S4, site_gen, site_irrep  )

site_irrep = { 0: Id }
ebr( S4, site_gen, site_irrep  )
'''

# point group C2 

'''
C2 = getPG(3)
for op in C2.op :
	print( op.name )

site_gen = [0] # Id
Id = sigma0
site_irrep = { 0: Id }
ebr( C2, site_gen, site_irrep )
'''

# point group \bar{3} ( C3i )

'''
C3i = getPG(17)
for op in C3i.op :
	print( op.name )

#exit()

site_gen = [0] # Id
Id = sigma0
site_irrep = { 0: Id }
ebr( C3i, site_gen, site_irrep )
'''

# point group C2h

'''
C2h = getPG(5)
for op in C2h.op :
	print( op.name )

#exit()

#site_gen = [1] # C2
site_gen = [3] # M
Id = sigma0
#C2 = exp_su2( 1, [0,0,1], 2, 1 )
M = exp_su2( 1, [0,0,1], 2, 1 )
site_irrep = { 0: Id, 3: M  }
ebr( C2h, site_gen, site_irrep )
'''

# point group 4mm

C4v = getPG(13)
for iop, op in enumerate(C4v.op):
	print( iop, op.name )
	print( np.round( op.su2, 1 ) )
	print( '\n' )

#print( C4v.get_coset( [0, 5] ) )
#print( C4v.mut[(3,6)] )
exit()

site_gen = [5] # Mx
Id = sigma0
Mx = complex(0,1) * sigma1
site_irrep = { 0: Id, 5: Mx }
ebr( C4v, site_gen, site_irrep )



