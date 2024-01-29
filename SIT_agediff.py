import os,asciidata
import numpy as np
import math
import decimal
from astLib import astCoords
from astLib import astCalc
from fractions import *

c=2.99792*10**5
PIG=2.118*10**11
Ms=4.68
MPC_km= 3.0857*10**19 #from MPC to meters

table= asciidata.open('SIT_z_age.dat')
fout = open('SIT_agediff.dat', 'w')
print>>fout, "#Index z12 z13  age12 age13 v12 v13  v_mean  d12 d13 d23 d_mean d_mean_2"

for i in range(table.nrows):
    Index = int(table[0][i])
    RA1   = float(table[1][i])
    DEC1  = float(table[2][i])
    RA2   = float(table[3][i])
    DEC2  = float(table[4][i])
    RA3   = float(table[5][i])
    DEC3  = float(table[6][i])
    z1    = float(table[7][i])
    z1_age= float(table[8][i])
    z2    = float(table[9][i])
    z2_age= float(table[10][i])
    z3    = float(table[11][i])
    z3_age= float(table[12][i])
    v1    = float(table[13][i])
    v2    = float(table[14][i])
    v3    = float(table[15][i])
    d1    = float(table[16][i])
    d2    = float(table[17][i])
    d3    = float(table[18][i])
    
    
    z12= abs(z1-z2)
    print "z12:", z12
    
    z13= abs(z1-z3)
    print "z13:", z13
    
    #z23= abs(z2-z3)
    #print "z23:", z23
    
    age12= abs(z1_age-z2_age)
    
    age13= abs(z1_age-z3_age)
    
    #age23= abs(z2_age-z3_age)
    
    v12= abs(v1-v2)
    
    v13= abs(v1-v3)
    
    #v23= (v2-v3)
    
    d12= abs(d1-d2)
    
    d13= abs(d1-d3)
    
    d23= abs(d2-d3)
 
    #v_mean= (v12+v13+v23)/3
    v_mean= (v12+v13)/2
    d_mean_2=(d12+d13)/2
    d_mean= (d12+d13+d23)/3
    
    print>>fout, Index, z12,z13, age12, age13, v12, v13, v_mean, d12, d13,d23, d_mean, d_mean_2