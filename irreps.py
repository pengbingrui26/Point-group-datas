# =================================================

import numpy as np
#from numpy import * 
from all_pg_irreps import all_pg_irreps
from PG import symop, Pointgroup
import pickle as pk
import os

# ======================================================

'''
for i in range( 32 ) :
	print( 'PG:', (i+1) )
	xx = all_pg_irreps[i]
	assert len( xx['dv_irreps'] ) == len( xx['sym_ops'] )  
	for j in range( len(xx['sym_ops']) ) :
		print( xx[ 'sym_ops' ][j][0] )
		for k in range( len(xx['dv_irrep_names']) ) :
			print( xx['dv_irrep_names'][k] )
			print( xx['dv_irreps'][j][k] )
		print( '\n' )
	print( '\n' )
'''

'''
for i in range( 32 ) :
	print( 'PG:', (i+1) )
	xx = all_pg_irreps[i]
	for op in xx[ 'sym_ops' ] :
		mat = op[1]
		find = False
		op_square = ''
		for opp in xx[ 'sym_ops' ] :
			mat_tmp = opp[1]
			if np.linalg.norm( np.dot( mat, mat ) - mat_tmp ) < 1e-6 :
				find = True 
				op_square = opp[0] 
		assert find 
		print( op[0], op_square )
	print( '\n' )
'''


# ==============================================================
# test rep_decompose

'''
from decompose_into_irreps import decompose_rep

for i in range( 31, 32 ) :
	print( 'PG:', (i+1) )
	xx = all_pg_irreps[i]
	assert len( xx['dv_irreps'] ) == len( xx['sym_ops'] )  
	reps = [  ]
	nrep = len( xx['dv_irrep_names'] )
	nop = len( xx['sym_ops'] )
	for i in range( nrep ) :
		rep = [ ]
		for j in range( nop ) :
			rep.append( xx['dv_irreps'][j][i] )
		reps.append( rep )
	
	#for rep in reps :
	#	print( rep )
	#	print( '\n' )
	print( len(reps), len(reps[0]) )
'''	

# ========================================================
# find all Kramers pairs 
# ========================================================

'''
for ipg in range( 32 ) :
	print( 'PG:', ipg+1 )
	pg = all_pg_irreps[ipg]
	op_square_dic = {  }
	nop = len( pg['sym_ops'] )
	for i in range( nop ) :
		opp = [ j for j in range(nop) if np.linalg.norm( \
			np.dot( pg['sym_ops'][i][1], pg['sym_ops'][i][1] ) \
			- pg['sym_ops'][j][1] ) < 1e-6 ]
		assert len(opp) == 1, ( len(opp), pg['sym_ops'][i][0] ) 
		opp = opp[0]
		op_square_dic[i] = opp
	#
	Kramers = [  ]
	nirreps = len( pg['dv_irrep_names'] )
	for ir in range( nirreps ) :
		irrep = [ pg['dv_irreps'][k][ir] for k in range(nop) ]
		irrep_name = pg['dv_irrep_names'][ir]
		#print( irrep_name )
		#print( irrep )
		assert len( irrep ) == nop		
		tot = 0
		for iop in range(nop) :
			mat = irrep[iop]
			mat_2 = irrep[op_square_dic[iop]]
			zgg = 0
			pl = np.linalg.norm( np.dot( mat, mat ) - mat_2 ) 
			mi = np.linalg.norm( np.dot( mat, mat ) + mat_2 ) 
			if pl < 1e-4 :
				zgg = 1
			elif pl > 1e-4 :
				zgg = -1
			assert zgg in [1,-1], ( ipg+1, irrep_name,  mat, mat_2 )
			assert zgg in [1,-1], ( pg['sym_ops'][iop][0], pg['sym_ops'][op_square_dic[iop]][0] )
			tr = np.trace( mat_2 ) 
			tot += zgg*tr
		tot = -tot/nop
		real = tot.real
		imag = tot.imag
		assert abs( real - round(real) ) < 1e-4
		assert abs( imag ) < 1e-4
		tot = int(round(real))
		#print( 'tot:', tot )
		assert tot in [ 1, -1, 0 ], tot 
		if tot == 1 :  # degenerace is unchanged
			pair = [ irrep_name ]
			if pair not in Kramers :
				Kramers.append( pair )
		elif tot == -1 : # two same irreps are paired
			pair = [ irrep_name, irrep_name ] 
			if pair not in Kramers :
				Kramers.append( pair )
		else :  # two conjugate irreps are paired 
			# find conjugate irrep
			irrep_conj = [ ma.conjugate() for ma in irrep ]
			assert len( irrep_conj ) == nop, len( irrep_conj )
			find_conj = False
			irrep_conj_name = ''
			for jr in range( nirreps ) :
				irrep_tmp = [ pg['dv_irreps'][k][jr] for k in range(nop) ]
				irrep_name_tmp = pg['dv_irrep_names'][jr]
				assert len( irrep_tmp ) == nop, len( irrep_tmp )
				if all( [ abs( np.trace(irrep_conj[i_tmp]) - np.trace(irrep_tmp[i_tmp]) ) \
						< 1e-4 for i_tmp in range(nop) ] ) :
					find_conj = True
					irrep_conj_name = irrep_name_tmp
					break
			assert find_conj, ( ipg+1, irrep_name )
			pair = sorted( [ irrep_name, irrep_conj_name ] ) 
			if pair not in Kramers :
				Kramers.append( pair )
		#print( '\n' )
	print( 'Kramers:', Kramers )	
	print( '\n' )
'''

# =====================================================
# Inverse (pair) 
# =====================================================

'''
for ipg in range( 32 ) :
	print( 'PG:', ipg+1 )
	pg = all_pg_irreps[ipg]
	op_square_dic = {  }
	nop = len( pg['sym_ops'] )
	for i in range( nop ) :
		opp = [ j for j in range(nop) if np.linalg.norm( \
			np.dot( pg['sym_ops'][i][1], pg['sym_ops'][i][1] ) \
			- pg['sym_ops'][j][1] ) < 1e-6 ]
		assert len(opp) == 1, ( len(opp), pg['sym_ops'][i][0] ) 
		opp = opp[0]
		op_square_dic[i] = opp
	#
	pairs = [  ]
	nirreps = len( pg['dv_irrep_names'] )
	for ir in range( nirreps ) :
		irrep = [ pg['dv_irreps'][k][ir] for k in range(nop) ]
		irrep_name = pg['dv_irrep_names'][ir]
		#print( irrep_name )
		#print( irrep )
		assert len( irrep ) == nop		
		find_inv = False
		irrep_inv_name = ''
		for jr in range( nirreps ) :
			irrep_tmp = [ pg['dv_irreps'][k][jr] for k in range(nop) ]
			irrep_name_tmp = pg['dv_irrep_names'][jr]
			assert len( irrep_tmp ) == nop, len( irrep_tmp )
			if all( [ abs( np.linalg.norm( irrep[i_tmp] + irrep_tmp[i_tmp] ) ) \
						< 1e-4 for i_tmp in range(nop) ] ) :
				find_inv = True
				irrep_inv_name = irrep_name_tmp
				break
		assert find_inv, ( ipg+1, irrep_name )
		pair = sorted( [ irrep_name, irrep_inv_name ] ) 
		if pair not in pairs :
			pairs.append( pair )
		#print( '\n' )
	print( 'pairs:', pairs )	
	print( '\n' )
'''

# =================================================================
# generate 32 PGs 
# =================================================================

'''
if not os.path.exists('PG') :
	os.mkdir('PG')		
	
for i in range( 32 ) :
	print( 'PG:', (i+1) )
	xx = all_pg_irreps[i]
	#print( xx['sym_ops'] )
	ops = [ symop( opp[1], opp[0] ) for opp in xx['sym_ops'] ]
	#for op in ops :
	#	print(op.mat)
	ipg = i+1
	pg = Pointgroup( ops, ipg )
	#print( pg.cls )
	#print( '\n' ) 
	name = './PG/'+str(ipg)+'.dat'		
	fd = open(name, 'wb')		
	pk.dump(pg,fd)		
	fd.close()
'''

# ========================================================================
# get point group 
# ========================================================================

def getPG( ipg ) :
	name = './PG/'+str(ipg)+'.dat'		
	fd = open(name, 'rb')		
	pg = pk.load(fd)		
	fd.close()
	return pg

# =====================================================================

T = getPG(28)
for op in T.op:
	print( op.name )
iC2x = 3
iC2z = 1
iC2y = 2 
iC3111p = 4

iC2x_inv = T.inv[iC2x]
iC2y_inv = T.inv[iC2y]
iC2z_inv = T.inv[iC2z]
#print( T.mut[ ( iC3111p, iC2z ) ] )
print( T.mut[ ( iC2x_inv, T.mut[ ( iC3111p, iC2z ) ] ) ] )
print( T.mut[ ( iC2y_inv, T.mut[ ( iC3111p, iC2x ) ] ) ] )
print( T.mut[ ( iC2z_inv, T.mut[ ( iC3111p, iC2y ) ] ) ] )
#subgroup = T.get_subgroup( [ iC2x,iC2y,iC2z ] ) 
subgroup = T.get_subgroup( [ iC3111p ] ) 
print( T.get_coset( subgroup ) )

'''
for i in range( 1, 33 ) :
	for op in getPG(i).op :
		print( op.name )
		print( op.su2 )
	print( '\n' )
'''

