import os,asciidata
import numpy as np
import math
import decimal
import matplotlib.pyplot as p
from astLib import astCoords
from astLib import astCalc
from scipy.stats.stats import pearsonr
from pylab import *

c=2.99792*10**5	  #km/sec
PIG=2.118*10**11  #9pi/2G 
Ms=4.68           #solar absolute magnitude
H0=70        	  # in km/S/Mpc
MPC_m= 3.0857*10**22 #from MPC to meters
Msun=1.9891*10**30 # solar mass in kg

 
a12_li=[]
a13_li=[]
a23_li=[]
rpsum_li=[]
Mvirl_li=[]


table2= asciidata.open('TGS_Mt2.dat')  #Index SIT RA DEC z kcorr mr vmean sigma sigmasq rs12 rs13 rs  rpmean LS Reff Mvirsun Mtrssun Mtrsl  Mtrpsun  Mtrpl

#fout= open('TGS_Mtrp.dat', 'w')
#print>>fout, "#Index rp Mtrp Mtrpl rp12 rp13 rp23"

f_out= open('TGS_rh_rs.dat','w')
print>>f_out,"#Index SIT RA DEC z sigma rh_rp rh_rs Mvir_rh R Mvir_rp Mvirl LS rp"

for j in range(table2.nrows):
    Index2 = int(table2[0][j])
    SIT   = int(table2[1][j])
    RA    = float(table2[2][j])
    DEC   = float(table2[3][j])
    z     = float(table2[4][j])
    kcorr = float(table2[5][j])
    mr    = float(table2[6][j])
    vmean = float(table2[7][j])
    sigma = float(table2[8][j])
    sigmasq= float(table2[9][j])
    rs12    = float(table2[10][j])
    rs13    = float(table2[11][j]) 
    rs23    = float(table2[12][j])
    rs      = float(table2[13][j])
    rpmean  = float(table2[14][j])
    LS      = float(table2[15][j])
    #Reff    = float(table2[15][j])
    #Mvirsun = float(table2[16][j])
    #Mtrssun = float(table2[17][j])
    #Mtrsl   = float(table2[18][j])
    #Mtrpsun = float(table2[19][j])
    #Mtrpl   = float(table2[20][j])
	
    
	
	    
    table = asciidata.open('SIT_zdiff.dat')


    ##Index sep_deg12 sep_arcmin12 sep_deg13 sep_arcmin13 sep_deg23 sep_arcmin23 z12 z13 z23

    for i in range(table.nrows):
	Index =    int(table[0][i])
	sep_deg12= float(table[1][i])
	sep_arcmin12= float(table[2][i])
	sep_deg13 = float(table[3][i])
	sep_arcmin13= float(table[4][i])
	sep_deg23 = float(table[5][i])
	sep_arcmin23 = float(table[6][i])
	z12 = float(table[7][i])
	z13= float(table[8][i])
	z23= float(table[9][i])	    
    
	
	if Index2==Index:
    
	    sep_rad12=math.radians(sep_deg12)
	    #print "seprad12:", sep_rad12
	    sep_rad13=math.radians(sep_deg13)
	    #print "seprad13:", sep_rad13
	    sep_rad23=math.radians(sep_deg23)
	    #print "seprad23:", sep_rad23
	    
	    y12=math.sin(sep_rad12/2.0)
	    #print "y12:", y12
	    y13=math.sin(sep_rad13/2.0)
	    #print "y13:", y13
	    y23=math.sin(sep_rad23/2.0)
	    #print "y23:", y23
	    
	    rp12=(2*vmean*y12)/H0
	    #print "rp12:", rp12
	    rp13=(2*vmean*y13)/H0
	    #print "rp13:", rp13
	    rp23=(2*vmean*y23)/H0
	    #print "rp23:", rp23
	    
	    rp=(rp12+rp13+rp23)/3
	    print "rp:", rp
	    
	    a1=math.pow(rp12,-1)
	    print "a1:", a1
	    a2=math.pow(rp13,-1)
	    print "a2:", a2
	    a3=math.pow(rp23,-1)
	    print "a3:", a3
	    
	    a=(a1+a2+a3)
	    print "a:", a
	    
	    x1=math.pow(rs12,-1)
	    x2=math.pow(rs13,-1)
	    x3=math.pow(rs23,-1)
	    x=(x1+x2+x3)/3
	    
	    
	    #rpsum_li.append(rpsum)
	    #print "rpsum_li:", rpsum_li
	    
	    #rp=sum(rpsum_li)/len(rpsum_li)
	    #print "rp:", rp
	    
	    
	    
	    # to get the total mass w.r.t. rp
	    
	    Mtrp=(PIG*(rp*MPC_m)*sigmasq*10**6)/Msun
	    #print "Mtrp:", Mtrp
	    
	    # to get M/L
	    
	    
	    Mtrpl= Mtrp/(LS)
	    #print "Mtrpl:", Mtrpl
	    
	    ##mean harmonic separation rh=[1/3(sum1/rp)]-1
		    
	    #a= math.pow(rp,-1)
	    b= a/3
	    print "b:", b
	    rh_rp= math.pow(b,-1)
	    #print "rh:", rh
	    
	    rh_rs=math.pow(x,-1)
	    
	    
	    Mvir_rh=(PIG*(rh_rp*MPC_m)*sigmasq*10**6)/Msun
	    print "Mvirrh:", Mvir_rh
	    
	    if rp12>rp13 and rp23:
		
		R=rp12/2
		
	    if rp13> rp12 and rp23:
		
		R=rp13/2
	    else:
		
		R=rp23/2
		
	    Mvir_rp=(PIG*(R*MPC_m)*sigmasq*10**6)/Msun
	    print "Mvir_rp:", Mvir_rp
	    
	    
	    Mvirl= Mvir_rh/LS
	    print "Mvr/l:", Mvirl
	    
	    print>>f_out,Index2,SIT, RA, DEC, z,sigma,rh_rp, rh_rs, Mvir_rh, R, Mvir_rp, Mvirl, LS, rp 
   
	    #print>>fout,Index, rp, Mtrp, Mtrpl,rp12, rp13, rp23
	    
	    Mvirl_li.append( Mvirl)
	    Mvirl_array=np.asarray(Mvirl_li)
	    
mMl= mean(Mvirl_array)
print "meanMvirl:", mMl

mdMl= median(Mvirl_array)
print "medianMvirl:", mdMl

stdMl= np.std(Mvirl_li, ddof=1)
print "stdMl:", stdMl
	    
	    
	    
   

f_out.close()
