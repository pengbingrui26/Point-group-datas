import numpy as np
import re
import pickle as pk
import os		
from all_pg_irreps import all_pg_irreps

# =============================================================
sigma0 = np.array( [ [ 1,0 ], [0,1] ] ).astype( complex )
sigma1 = np.array( [ [ 0,1 ], [1,0] ] ).astype( complex )
sigma2 = np.array( [ [ 0,-complex(0,1) ], [complex(0,1),0] ] \
						).astype( complex )
sigma3 = np.array( [ [ 1,0 ], [0,-1] ] ).astype( complex )
sigma = [ sigma0, sigma1, sigma2, sigma3 ]
# ==============================================================

def is_number(s):		
	try:		
		float(s)		
		return True		
	except ValueError:		
		pass		
		
	try:		
		import unicodedata		
		unicodedata.numeric(s)		
		return True		
	except (TypeError, ValueError):		
		pass		
	
	return False

class symop( object ) :
	def __init__( self, mat, name ) :
		self.mat = mat
		self.name = name
		# get SU(2) matrix
		self.su2 = ''
		#print( name[0], name[1], name[2] )
		assert len(name) == 3
		assert any( [ xx in name[0] for xx in \
						[ '1', 'm', '2', '3', '4', '6' ] ] ), name
		num = re.findall( r'\d+\.?\d*', name[0] ) 
		assert len(num) in [0,1], num
		pm = ''
		if '+' in name[1] :
			pm = 1
		elif '-' in name[1] :
			pm = -1
		else :
			assert '1' in name[0] or '2' in name[0] or 'm' in name[0]
			pm = 1
		assert pm in [1,-1], name
		#
		
		'''
		axis = ''
		if name[2] != '' :
			axis_tmp = [ ]
			for i, x in enumerate(name[2]) :
				if is_number( x ) :
					if i == 0 :
						axis_tmp.append( int(x) )
					else :  
						if name[2][i-1] == '-' :
							axis_tmp.append( -int(x) )
						else :
							assert is_number( name[2][i-1] ), name[2][i-1]
							axis_tmp.append( int(x) )
			assert len(axis_tmp) == 3, axis_tmp
			axis = axis_tmp
		if len(num) == 0 :
			assert 'm' in name[0], name
			assert axis != ''
			axis = axis/np.linalg.norm(axis)
			phi = np.pi/2
			axisdotsigma = axis[0] * sigma[1] + axis[1] * sigma[2] + axis[2] * sigma[3]
			self.su2 = np.cos(phi) * sigma0 + pm * complex(0,1) * np.sin(phi) * axisdotsigma
		elif len(num) == 1 :
			nfold = int(num[0])
			if nfold == 1 :
				self.su2 = sigma0
			else :
				assert axis != '', name
				axis = axis/np.linalg.norm(axis)
				phi = np.pi/nfold
				axisdotsigma = axis[0] * sigma[1] + axis[1] * sigma[2] + axis[2] * sigma[3]
				self.su2 = np.cos(phi) * sigma0 + pm * complex(0,1) * np.sin(phi) * axisdotsigma
		'''
		


class Pointgroup( object ) :
	def __init__( self, op, pgid ) :
		self.op = op
		self.pgid = pgid
		self.nop = len(self.op)
		self.opnames = [ xx.name for xx in self.op ]
		# identity element 
		id_mat = np.array( [ [1,0,0], [0,1,0], [0,0,1] ] )
		self.id = -1
		for i in range(self.nop) :
			if abs( np.linalg.norm( self.op[i].mat - id_mat ) ) < 1e-4 :
				self.id = i
				break
		assert self.id != -1
		# group multiplication table
		self.mut = {  }
		for i in range(self.nop) :
			for j in range(self.nop) :
				matij = np.dot( self.op[i].mat, self.op[j].mat )
				tmp = [ k for k in range(self.nop) if \
						abs( np.linalg.norm( matij-self.op[k].mat) ) < 1e-4 ]
				assert len(tmp) == 1, len(tmp)
				self.mut[(i,j)] = tmp[0]
		# inverse elements
		self.inv = { }
		for i in range(self.nop) : 
			inv = [ j for j in range(self.nop) if 
				abs( np.linalg.norm( np.dot( self.op[i].mat, \
											self.op[j].mat ) - id_mat ) ) < 1e-4 ]
			assert len(inv) == 1, len(inv)
			self.inv[i] = inv[0]
		# conjugate classes 
		self.cls = [ ]
		used = [ 'False' ] * self.nop
		for i in range(self.nop) :
			if used[i] == 'True' :
				continue
			cls = [ ]
			for j in range(self.nop) :
				k = self.mut[ ( self.inv[j], self.mut[(i,j)] ) ]
				if k not in cls:
					cls.append(k)
					used[k] = 'True'
			cls = sorted(cls)	
			if cls not in self.cls :
				self.cls.append(cls)

	def get_subgroup( self, subgen ) : # get the subgroup given a set of generators
		subops = subgen
		for i in subgen:
			for j in subgen:
				ij = self.mut[(i,j)]
				if ij not in subops:
					subops.append(ij)
		subops.sort()
		return subops

	def get_coset( self, subgroup ) :
		cosets = [ ]
		coset_reps = [ ]
		for iop in range(self.nop) :
			#tmp = [ self.mut[ ( subiop, iop ) ] for subiop in subgroup ]
			tmp = [ self.mut[ ( iop, subiop ) ] for subiop in subgroup ]
			tmp.sort()
			if tmp not in cosets :
				cosets.append(tmp)
				coset_reps.append(iop)
		coset_dict = { }
		for i in range(len(cosets)) :
			coset_dict[coset_reps[i]] = cosets[i]
		return coset_dict 

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
	#   print(op.mat)		
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

# ==========================================================================
# test PG 
# ==========================================================================

'''
S4 = getPG(10)
for op in S4.op :
	print( op.name )
	
print( S4.cls )
''' 
