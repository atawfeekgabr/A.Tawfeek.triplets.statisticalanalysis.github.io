
import os,asciidata
import numpy as np
from astLib import astCoords
from astLib import astCalc

table1 = asciidata.open('TGS.dat')

table2= asciidata.open('tab1.dat')

c=2.99792*10**5
PIG=2.118*10**11
Ms=4.68

#for i in range(table1.nrows):
    #Index1 = table1[0][i]
    #SIT1   = table1[1][i]
    #RA1    = float(table1[2][i])
    #DEC1   = float(table1[3][i])
    #z1     = float(table1[4][i])
    #kcorr1 = float(table1[5][i])
    #mr1    = float(table1[6][i])
    
    #if int(Index1) == 1:
	
	#print Index1, RA1, DEC1    

#v_li = []

for j in range(table2.nrows):
    Index2 = table2[0][j]
    SIT2   = table2[1][j]
    RA2    = float(table2[2][j])
    DEC2   = float(table2[3][j])
    z2     = float(table2[4][j])
    kcorr2 = float(table2[5][j])
    mr2    = float(table2[6][j])
    v      = float(table2[7][j])
    D      = float(table2[8][j])
    LD     = float(table2[9][j])
    Mr     = float(table2[10][j])
    Mrk    = float(table2[11][j]) 
    L      = float(table2[12][j])
    
    #x= table2[7][0]+table2[7][1]+table2[7][2]
    #print x
    
    v_li = [table2[7][0],table2[7][1],table2[7][2]]
    print v_li
    
    x=sum(v_li)
    print x
    
    y= sum(v_li)/len(v_li)
    print y
    
    z= y*y
    print z
    
    v2_li = [table2[7][0]*table2[7][0],table2[7][1]*table2[7][1],table2[7][2]*table2[7][2]]
    print v2_li
    
    sv= sum(v2_li)/ len(v2_li)
    print sv
    
    sigmasq= sv-z
    print sigmasq
    
    sigma= np.sqrt(sigmasq)
    print sigma
    
    #v_li.append(v)
    
	    
	    #print np.size(v)
	    #print v, np.mean(v)
	
	    #sigma = float (np.sqrt(np.sum((v - np.mean(v))**2) / np.size(v)))
	
	    #print("Velocity dispersion: {0:.2f}".format(sigma))
	
	#v_li.append(v)
	    ##v**2_li.append(v**2)
	    
	#print len(v_li)
	
	    #if RA1==RA2:	
	
	    #sep_deg12    =  astCoords.calcAngSepDeg(RA1, DEC1, RA2, DEC2)
	    #sep_arcmin12 =  float(sep_deg12*60.0) 
	    #sep_arcsec12 =  float(sep_deg12*3600)
	    
	    #print  "sep_arcmin12:", sep_arcmin12
	    #print Index1, Index2, SIT2,  sep_arcmin12
	    #print sep_arcmin13
	
	#D = v/70.0
	
	#LD = astCalc.dl(z2)

	#Mr = astCalc.absMag(mr2,LD)
	
	#Mrk = Mr-kcorr2
	
	#X = Ms-Mrk
	
	#Y = 0.4*X
	
	#L  = 10**Y
	
			

