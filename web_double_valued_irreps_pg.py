# =================================================================
# get all double-valued irreps from Bilbao
# =================================================================
from selenium import webdriver
import time
import pickle as pk

wb = webdriver.Chrome()

C4page = 'https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_out.pl?tipogrupo=dbg&pointspace=point&num=75&super=9&symbol=4'

main_page = 'https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations_point.pl?tipogrupo=dbg'

wb.get( main_page )
#wb.get( C4page )

# =======================================================
# whole source page 
#

#page = wb.page_source 
#print( page )

# =======================================================

selections = wb.find_elements_by_xpath( " //*[@bgcolor='#f5f5f5'][@align='center'] " )

tables = { }
#tables = [ ]

i = 0
#while i < 32 :
while i < 2 :
#while i < len( selections ) :
    selection = wb.find_elements_by_xpath( " //*[@bgcolor='#f5f5f5'][@align='center'] " )[i]
    selection.click()
    time.sleep(2)
    #table = wb.find_element_by_xpath( " //*[@border='5']/tbody " )
    #tables[ i ] = table.get_attribute( 'innerHTML' )
    datas = [ ] 
    #trlist = table.find_elements_by_tag_name( 'tr' )
    #trlist = table.find_elements_by_css_selector( ' #tr:nth-child(1) ' )
    trlist = wb.find_elements_by_xpath( " //*[@border='5']/tbody/tr " )
    #print( 'len_trlist:', len(trlist) )
    for row in trlist :
        #print( row.get_attribute('innerHTML') )
        data_rows = [ ]
        #tdlist = row.find_elements_by_tag_name( 'td' )
        #tdlist = row.find_elements_by_css_selector( ' td:nth-child(1) ' )
        tdlist = row.find_elements_by_xpath( 'td' )
        #print( 'len_tdlist', len(tdlist) )
        for col in tdlist :
            #print( col.get_attribute('innerHTML') )
            data_rows.append( col.get_attribute('innerHTML') ) 
            print( col.get_attribute('innerHTML') ) 
            print( '\n' )
        datas.append( data_rows )
        #print( data_rows )
        #print( '\n' )
    tables[i+1] = datas
    #tables.append( datas )
    wb.back()
    i += 1

#for k in tables :
#    print( tables[k] )
#    print( '\n' )

'''
f = open( './datas.txt', 'wb' )
pk.dump( tables, f )
f.close()
'''



