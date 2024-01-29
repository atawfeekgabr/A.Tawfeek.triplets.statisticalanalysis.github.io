import os,asciidata
import numpy as np
from astLib import astCoords
from astLib import astCalc

table1 = asciidata.open('TGS.dat')

table2= asciidata.open('SIT_NEW.dat')
#fout = open('TGS2.dat', 'w')
 #print>>fout,'#Index, SIT, RA, DEC, z, kcorr, mr'

for i in range(table1.nrows):
    Index1 = table1[0][i]
    SIT1   = table1[1][i]
    RA1    = table1[2][i]
    DEC1   = table1[3][i]
    z1     = float(table1[4][i])
    kcorr1 = float(table1[5][i])
    mr1    = float(table1[6][i])
    
    #fil = open("'tab%s.dat'", 'w' %Index1)
    fil = open('tab{0}.dat'.format(Index1), 'w')

    print>>fil, '#Index SIT RA DEC z kcorr mr v D LD Mr Mrk L'
    
    print i+1,   Index1, SIT1, RA1, DEC1, z1, kcorr1, mr1
    
    for j in range(table2.nrows):
	Index2 = table2[0][j]
	SIT2   = table2[1][j]
	RA2    = table2[2][j]
	DEC2   = table2[3][j]
	z2     = float(table2[4][j])
	kcorr2 = float(table2[5][j])
	mr2    = float(table2[6][j])
	v2     = float(table2[7][j])
	D2     = float(table2[8][j])
	LD2    = float(table2[9][j])
	Mr2    = float(table2[10][j])
	Mrk2   = float(table2[11][j])
	L2     = float(table2[12][j])
	
	if Index2==Index1 :
	    
	    #print Index1,  Index2
	    print>>fil, Index2, SIT2, RA2, DEC2, z2, kcorr2, mr2, v2, D2, LD2, Mr2, Mrk2, L2
      
      
    fil.close()
    
    #os.system('mv %s.dat tables_TS/' %Index1)
    #os.system('mv tab{0}.dat tables_TS/'.format(Index1))