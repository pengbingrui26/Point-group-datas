# ====================================================================
# read data.txt to extract double-valued irreps
# ====================================================================

import pickle as pk
import re
import numpy as np

f = open( './datas.txt', 'rb' )
datas = pk.load( f )
f.close()

# decide whether a str is a number ==================================

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

# datas for all 32 point groups ======================================================
# a dict for each point group: {  double-valued irrep names: xxx, symm_op names: xxx, \
#				            double-values irreps : [ [xx], [xx], ...., [xx] ]  } 
# ===================================================================================

all_32_datas = [  ]
for i in range( 32 ) :
	all_32_datas.append( { } )  


# ===============================================================================
# try dealing with double-valued irrep names 
# ===============================================================================

for ipg in range( 1, 33 ) :
	data = datas[ipg]
	assert len(data)%2 == 1, len(data)
	names = [ ]
	for dd in data[0] :
		if 'overline' in dd :
			assert 'GM' in dd, dd	
			index = re.findall( r'<sub>(.+?)</sub>', dd )	
			prop = re.findall( r'</sub>\((.+?)\)</center>', dd )		
			index, prop = ''.join( index ), ''.join( prop )
			#names.append( 'GM' + ' overline' + ' ' + index + ' ' + prop )
			names.append( 'GM' + index )
	all_32_datas[ipg-1]['dv_irrep_names'] = names
	#print( '\n' )

'''
for da in all_32_datas :
	print( da )
'''

# try dealing with sym_op names and SO(3) mats ===================================

for ipg in range( 1, 33 ) :

	data = datas[ipg]
	assert len(data)%2 == 1, len(data)
	sym_ops = [  ]

	for k in range( 1, int((len(data)+1)/2) ) :
		strs = data[k][3]
		mat_strs = data[k][1]
		
		## deal with sym_op name :
		index_tmp = re.findall( r'<center>(.+?)<', strs )
		index = ''
		if re.search( 'overline', index_tmp[0] ) :
			num = re.findall( r'\d+\.?\d*', index_tmp[0] ) 
			assert len(num) == 1, len(num)
			index = '\\bar{' + num[0] + '}'
		else :
			index = index_tmp[0]

		pm = ''	
		#if re.search( '>+<', strs ) :
		if '<sup>+</sup>' in strs :
			pm = '+'
		#elif re.search( '<sup>-', strs ) :
		elif '<sup>-</sup>' in strs  :
			pm = '-'

		axis = ''
		axis_tmp = re.findall( r'<sub>(.+?)</sub>', strs )
		assert len(axis_tmp) in [ 0,1 ], axis_tmp
		if len(axis_tmp) == 1 :
			axis_tmp = axis_tmp[0]
			if 'overline' not in axis_tmp :
				nums = re.findall( r'\d+\.?\d*', axis_tmp ) 
				nums = ''.join( nums )
				axis = nums
			elif 'overline' in axis_tmp :
				assert sum( 1 for xx in axis_tmp if xx == 's' ) in [1,2,3], axis_tmp
				num_index = [ ]
				overline_index = [ ]
				index_dic = {  } 
				for k in range( len(axis_tmp) ) :
					st = axis_tmp[k]
					if is_number( st ) :
						num_index.append( k )
						index_dic[k] = st
					elif st == 's' : 
						index_dic[k] = '-'
						overline_index.append( k )
				assert len( num_index ) == 3, axis_tmp
				assert len( overline_index ) in [1,2,3], axis_tmp
				for ii in sorted( index_dic ) :
					axis += index_dic[ii] 

		## deal with sym_op SO(3) mat :
		mat_elemts = re.findall( r'<td align=\"right\">(.+?)</td>', mat_strs )
		assert len(mat_elemts) == 9, ( len( mat_elemts ), mat_strs ) 
		mat = [ ] 
		for elem in mat_elemts :
			nums = re.findall( r'\d+\.?\d*', elem )
			assert len(nums) == 1, ( len(nums), nums )
			num = nums[0] 
			if '-' in elem :
				mat.append( -int(num) )
			else : 
				mat.append( int(num) )
		len_mat = len(mat)
		assert len_mat == 9, mat
		mat = np.array( mat ).reshape( 3, 3 )
		
		sym_ops.append( [ ( index, pm, axis ), mat ] )
		
	all_32_datas[ipg-1]['sym_ops'] = sym_ops
	

# try dealing with irreps =========================================

def deal_exp( strs ) :
	ia = re.findall( r'e<sup>(.+?)π', strs ) 
	b = re.findall( r'π/(.+?)</sup>', strs )
	assert len(ia) == 1, ia
	assert len(b) == 1, b
	ia, b = ''.join( ia ), ''.join( b )
	nums = re.findall( r'\d+\.?\d*', ia )
	minus = ''
	minus_value = 1
	if re.search( '-', ia ) != None :
		minus = '-'
		minus_value = -1
	#print( 'minus:', minus )
	a = ''
	a_value = 1
	if len(nums) != 0 : 
		assert len(nums) == 1
		a = nums[0]
		a_value = int(a)
	b = int(b)
	#return a, b
	#return np.exp( minus * complex(0,1) * np.pi * a / b )
	#print( type(minus_value), type(a), type(b) )  
	return 'exp' + minus + 'iπ' + str(a) + '/' + str(b) ,  \
			np.exp( minus_value * complex(0,1) * np.pi * a_value / b )


def deal_plus_i( strs ) :
	nums = re.findall( r'\d+\.?\d*', strs ) 
	assert len(nums) == 1, len(nums) 
	return int( nums[0] )

def deal_minus_i( strs ) :
	nums = re.findall( r'\d+\.?\d*', strs ) 
	assert len(nums) == 1, len(nums) 
	return int( nums[0] )

def deal_sqrt( strs ) :
	#nums = re.findall( r'\d+\.?\d*', strs ) 
	nums = re.findall( r'text-decoration:overline;\">(.+?)</span>', strs )
	assert len(nums) == 1, ( len(nums), strs )
	#return int( nums[0] )
	return nums[0], int(nums[0])

def deal_divide( strs ) :
	nums = re.findall( r'\d+\.?\d*', strs ) 
	assert len(nums) == 1, ( len(nums), strs ) 
	#return int( nums[0] )
	return nums[0], int( nums[0] )

def deal_elem( strs_init ) :
	#print( 'strs:', strs_init )
	strs_1 = strs_init.lstrip()
	strs = strs_1.rstrip()
	
	ia = re.findall( r'e<sup>(.+?)π', strs ) 
	b = re.findall( r'π/(.+?)</sup>', strs )
	#plus_i = re.findall( r'i<span', strs )
	#minus_i = re.findall( r'-i<span', strs )
	#plus_1 = re.findall( r'<span', strs )
	#minus_1 = re.findall( r'-<span', strs )

	plus_i = strs.startswith( r'i<span' )
	minus_i = strs.startswith( r'-i<span' )
	plus_1 = strs.startswith( r'<span' )
	minus_1 = strs.startswith( r'-<span' )
	
	sqrt = re.findall( r'nowrap">√<span', strs )
	divide = re.findall( r'>(.+?)</span></span>/', strs )
	
	assert ( len(ia) == len(b) ), strs
	len_exp = len(ia)
	assert len_exp in [0,1], strs
	
	assert all( len(xx) in [0,1] for xx in [ sqrt, divide ] ), strs
	
	result = [  ]
	result_value = 1
	
	if len( strs ) <= 10 :
		if 'i' in strs :
			assert len( re.findall( r'\d+\.?\d*', strs ) ) == 0, strs
			if '-' in strs :
				result.append( '-i' )
				result_value *= complex(0,1) 
			else :
				result.append( 'i' )
				result_value *= -complex(0,1) 
		else :
			assert len( re.findall( r'\d+\.?\d*', strs ) ) == 1, strs
			if '-' in strs :
				result.append(  '-' + re.findall( r'\d+\.?\d*', strs )[0] )
				result_value *= -1 * int( re.findall( r'\d+\.?\d*', strs )[0] )
			else :
				result.append(  re.findall( r'\d+\.?\d*', strs )[0] )
				result_value *= int( re.findall( r'\d+\.?\d*', strs )[0] )
	
	else : 
		if len_exp == 1 :
			#assert all( len(xx) == 0 for xx in [ plus_i, minus_i, plus_1, minus_1 ] ), strs
			assert all( xx == False for xx in [ plus_i, minus_i, plus_1, minus_1 ] ), strs
			result.append( deal_exp( strs )[0] )
			result_value *= deal_exp( strs )[1]
			if len( sqrt ) == 1 :
				result.append( '√' + deal_sqrt( strs )[0] )
				result_value *= np.sqrt( deal_sqrt( strs )[1] )
			if len( divide ) == 1 :
				result.append( '/' + strs[-1] )
				result_value /= float( int( strs[-1] ) )
		
		else :
			tot = [ xx for xx in [ plus_i, minus_i, plus_1, minus_1 ] if xx == True ]
			assert len(tot) in [0,1], strs
			if len(tot) == 1 :
				assert len(sqrt) == 1 and len(divide) == 1, strs
				if tot[0] == plus_i :
					result.append( 'i' )
					result_value *= complex(0,1)
				elif tot[0] == minus_i :
					result.append( '-i' )
					result_value *= -complex(0,1)
				#elif tot[0] == plus_1 :
				#    result.append( '' )
				elif tot[0] == minus_1 :
					result.append( '-' )
					result_value *= -1
				result.append(  '√' + deal_sqrt( strs )[0] )
				result_value *= np.sqrt( deal_sqrt( strs )[1] )
				result.append(  '/' + deal_divide( divide[0] )[0] )
				result_value /= float( deal_divide( divide[0] )[1] )
			else :
				assert len( sqrt ) == 0 and len( divide ) == 0, strs
 	
	return ''.join( result ), result_value


# =================================================================================

for ipg in range( 1, 33 ) :
	data = datas[ipg]
	assert len(data)%2 == 1, len(data)

	## find the begining of double-valued irreps 
	begin_dv_irrep = 0
	find = False
	#while find == False :
	for i in range( len(data[0]) ):
		dd = data[0][i]
		if 'overline' in dd :
			begin_dv_irrep = i
			assert 'GM' in dd, dd	
			find = True
			break
	#print( begin_dv_irrep, len(data[0]), len(data[1]) )

	dv_irreps = [ ]	
	dv_irreps_str = [ ]
	for k in range( 1, int((len(data)+1)/2) ) :
		dv_irreps_symmop = [ ]
		dv_irreps_symmop_str = []
		#print( 'begin_dv_irrep, len_data_k :', begin_dv_irrep, len(data[k]) )
		#print( 'data[k][-1]:', data[k][-1] )
		for m in range( begin_dv_irrep+1, len(data[k]) ) : 
			strs = data[k][m]
			tbs = re.findall( r'<tbody>(.+?)</tbody>', strs )
			assert len( tbs ) == 1, strs
			tbs = ''.join( tbs )
			bb = re.findall( r'align=\"left\">(.+?)</td>', tbs ) + \
					re.findall( r'align=\"right\">(.+?)</td>', tbs ) 
			#if k == 1 :
			#    len_bb = len(bb)
			#else :
			#    assert len_bb == len(bb), ( len_bb, len(bb) )
			#    len_bb = len(bb)
			#print( 'len(bb):', len(bb) )
			#if len(bb) != 0 :
			#if len(bb) == 4 :
				#print( len(bb) )
				#print( bb )
				#assert len(bb) in [1, 4, 16], len(bb)
			mat_elemts_str = [  ]
			mat_elemts = [  ]
			for bbb in bb :
				assert len( deal_elem(bbb) ) <= 20
                #print( '\n' )
				mat_elemts_str.append( deal_elem(bbb)[0] )
				mat_elemts.append( deal_elem(bbb)[1] )
				#if 'e<sup>' in bbb and 'π' in bbb :
					#print( bbb )
					#print( deal_exp( bbb ) )
					#mat_elemts.append( deal_exp(bbb) )
				#else :
					#mat_elemts.append( bbb )
			assert len(mat_elemts_str) in [ 1, 4, 9, 16 ], len(mat_elemts_str)
			num_elemts = len( mat_elemts )
			assert num_elemts in [ 1, 4, 9, 16 ], num_elemts
			mat_elemts = np.array( mat_elemts).astype( complex )  
			mat_elemts = mat_elemts.reshape( int( round(np.sqrt(num_elemts) ) ), -1 )
			dv_irreps_symmop_str.append( mat_elemts_str )
			dv_irreps_symmop.append( mat_elemts )
			#print( '\n' )
		dv_irreps_str.append( dv_irreps_symmop_str )
		dv_irreps.append( dv_irreps_symmop )
		#print( '\n' )
	#print( '\n' )
	#all_32_datas[ipg-1]['dv_irreps_str'] = dv_irreps_str
	all_32_datas[ipg-1]['dv_irreps'] = dv_irreps

print( 'import numpy as np' )
print( 'all_pg_irreps = ', all_32_datas )


'''
#for da in all_32_datas:
for i in range( 30, 31 ) :
	da = all_32_datas[i]
	#print( len( da['dv_irreps'] ), len( da['dv_irreps'][0] ) )
	assert len( da['dv_irreps'] ) == len( da['symmop_names'] )
	assert len( da['dv_irreps'][0] ) == len( da['dv_irrep_names'] )
	#for aa in da['dv_irreps'] :
	#	for aaa in aa:
	#		print(aaa)
	#for aa in da['dv_irreps'][-1] :
	#	print(aa)
	for j in range( len(da['dv_irreps'][-1]) ) :
		print( da['dv_irreps_str'][-1][j] )
		print( da['dv_irreps'][-1][j] )
	print( '\n' )
'''

# =================================================================

class irrep_name( object ) :
    def __init__( self, index, prop, sd ) :  
        # i.e., GM6(0) has index 6, and prop 'real,' 
        # where '1' for 'real', '-1' for 'pseudoreal', and '0' for 'complex'
        # sd means for spinless or spin-1/2
        self.index = index
        self.prop = prop
        self.sd = sd
        if self.sd == 's' :
            self.name_str = 'GM' + str(index) + '(' + str(prop) + ')' + 's' 
        elif self.sd == 'd' :
            self.name_str = 'GM' + str(index) + '(' + str(prop) + ')' + 'd' 

class sym_name( object ) :
    def __init__( self, index, pm, axis ) :  
        # index can be 1~6 and 'm'
        # pm for '+' or '-'
        self.nfold = nfold
        self.pm = pm
        self.axis = axis



