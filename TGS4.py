import os,asciidata
import numpy as np
import math
import decimal
from astLib import astCoords
from astLib import astCalc
from fractions import *


c=2.99792*10**5	  #km/sec
PIG=2.118*10**11  #9pi/2G 
Ms=4.68           #solar absolute magnitude
H0=70        	  # in km/S/Mpc
Msun=1.9891*10**30 # solar mass in kg

table1 = asciidata.open('TGS_sigma.dat')
fout= open('TGS_Mt2.dat', 'w')
print>>fout, "#Index SIT RA DEC z kcorr mr vmean sigma sigmasq rs12 rs13 rs23 rs  rpmean LS Reff Mtrssun Mtrsl  Mtrpsun  Mtrpl Mvir_rs"


for i in range(table1.nrows):
    Index1 = int(table1[0][i])
    SIT1   = int(table1[1][i])
    RA1    = float(table1[2][i])
    DEC1   = float(table1[3][i])
    z1     = float(table1[4][i])
    kcorr1 = float(table1[5][i])
    mr1    = float(table1[6][i])
    sigma  = float(table1[7][i])
    sigmasq= float(table1[8][i])
  

    table2= asciidata.open('table%s.dat' %Index1)
    f_out = open('table_new%s.dat' %Index1, 'w')
    print>>f_out, "#Index SIT RA DEC z kcorr mr v D LD Mr Mrk L sep_deg  sep_arcmin sep_arcsec rp "


    sep_deg_li =[]
    sep_arcmin_li =[]
    sep_rad_li = []
    a_li = []
    rp_li = []
    l_li=[]


    for j in range(table2.nrows):
	Index2 	   = int(table2[0][j])
	SIT2       = int(table2[1][j])
	RA2        = float(table2[2][j])
	DEC2       = float(table2[3][j])
	z2         = float(table2[4][j])
	kcorr2     = float(table2[5][j])
	mr2        = float(table2[6][j])
	v          = float(table2[7][j])
	D          = float(table2[8][j])
	LD         = float(table2[9][j])
	Mr         = float(table2[10][j])
	Mrk        = float(table2[11][j]) 
	L          = float(table2[12][j])
	sep_deg    = float(table2[13][j])
	sep_arcmin = float(table2[14][j])
	sep_arcsec = float(table2[15][j])
	
	
	sep_deg_li.append(sep_deg)
	sep_arcmin_li.append(sep_arcmin)
	
	l_li.append(L)
	LS= sum(l_li)
	
	v_li = [table2[7][0],table2[7][1],table2[7][2]]
	#print v_li

	
	
	if sep_deg > 0:
	    
	    
	    sep_rad=math.radians(sep_deg)
	    #print "sep_rad:", sep_rad
	
	    #sep_rad_li.append(sep_rad)
	    #print "sep_rad:", sep_rad_li
	    
	    x= math.cos(sep_rad)
	    #print x
	    
	    y= math.sin(sep_rad/2.0)
	    #print y
	    
	    z= sum(v_li)/len(v_li)
	    #print z
	    
	    rpH= 2*z*y
	    rp= rpH/H0
	    #print "rp:",rp
	    
	    rp_li.append(rp)
	    rpmean= sum(rp_li)/len(rp_li)
	    #print "rpmean:", rpmean
	    
	    rpmean_m = rpmean*3.0857*10**22 #from MPC to m
	    
	    Mtrp= PIG*rpmean_m*(sigmasq*10**6)
	    Mtrpsun= Mtrp/Msun
	    #print "Mvirrpsun:", Mtrpsun
	    
	    # luminosity l=dex[0.4(Mrsun - Mr)], dex is 10**, Mrsun is the abs. mag. of sun in r-band, Mr is the abs.mag.of g.in r-band.
	    #L in table was L/10**10, so I have to multiply all L by 10**10
		
	    Mtrpl= Mtrpsun/(LS)
	    #print "Mtrp/L:", Mtrpl
	    
	    
	    ##mean harmonic separation rh=[1/3(sum1/rp)]-1
	    
	    a= math.pow(rp,-1)
	    #print "a:", a
	    
	    a_li.append(a)
	    
	    b= sum(a_li)/len(a_li)
	    #print "b:", b
	    
	    rh= math.pow(b,-1)
	    print "rh:", rh
	    
	    
	    v12= table2[7][0]*table2[7][1]
	    #print v12
	    
	    v13= table2[7][0]*table2[7][2]
	    #print v13
	    
	    v23=table2[7][1]*table2[7][2]
	    

	    
	    rsH12= np.sqrt(((table2[7][0]**2) + (table2[7][1]**2)) - (2*v12*x))
	    rs12= rsH12/H0
	    #print "rs12:", rs12
	    
	    rsH13= np.sqrt(((table2[7][0]**2) + (table2[7][2]**2)) - (2*v13*x))
	    rs13= rsH13/H0
	    #print rs13
	    
	    rsH23= np.sqrt(((table2[7][1]**2) + (table2[7][2]**2)) - (2*v23*x))
	    rs23= rsH23/H0
	    print "rs23:", rs23
	    
	    rs=(rs12+rs13+rs23)/3    #to calculate <rs>
	    print "rs:",rs
	    
	    #if rs > 1.0:
		
		#rsF=rs
		#print "rsF:", rsF
		
	    #else:
		#rsT=rs
		#print "rsT:", rsT
	    
	    rs_m=rs*3.0857*10**22
	    
	    Mtrs=PIG*rs_m*(sigmasq*10**6)
	    Mtrssun=Mtrs/Msun
	    #print "Mtrs:",Mtrssun
	    
	    Mtrsl=Mtrssun/(LS)
	    #print "Mtrsl:", Mtrsl
	    
	
		
	    if rs12>rs13 and rs23:
	    
		Reff= rs12/2
	    
	    if rs23>rs12 and rs13:
	       
	       Reff= rs23/2
	       
	    else:
		Reff = rs13/2
		#print "Reffective:", Reff
		
	        #from MPC to meters
	        Reffmeters = Reff*3.0857*10**22
	        #print "Reffmeters:", Reffmeters
	        
	        sigma_meters= sigmasq*10**6
	        
		Mvir=PIG*Reffmeters*sigma_meters
		#print Mvir
	
		Mvirsun= Mvir/Msun
		print "Mvir/Msun:", Mvirsun

	
	
	else :
	    rp=0	    
	    
	    
		
	    print j+1,    Index2, SIT2, RA2, DEC2, z2, kcorr2, mr2, v, D, LD, Mr, Mrk, L, sep_deg, sep_arcmin, sep_arcsec, rp
	    print>>f_out, Index2, SIT2, RA2, DEC2, z2, kcorr2, mr2, v, D, LD, Mr, Mrk, L, sep_deg, sep_arcmin, sep_arcsec, rp
		
    print i+1,  Index1, SIT1, RA1, DEC1, z1, kcorr1, mr1, z, sigma, sigmasq, rs12, rs13, rs23, rs,rpmean,LS, Reff, Mtrssun, Mtrsl, Mtrpsun,  Mtrpl, Mvirsun 
    print>>fout,Index1, SIT1, RA1, DEC1, z1, kcorr1, mr1, z, sigma, sigmasq, rs12, rs13, rs23, rs,rpmean,LS, Reff, Mtrssun, Mtrsl, Mtrpsun,  Mtrpl, Mvirsun

f_out.close()

fout.close()