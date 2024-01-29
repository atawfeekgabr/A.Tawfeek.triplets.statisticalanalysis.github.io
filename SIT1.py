#to read table
import os,asciidata
from astLib import astCoords
from astLib import astCalc
from math import *
import numpy as np 

#print pi

table = asciidata.open('SIT_input5_new.dat')

c=2.99792*10**5
#PIG=9pi/2G (constant for virial mass equation)
PIG=2.118*10**11
# solar absolute mag in r-bad
Ms=4.68
#print c

fout = open('SIT_NEW.dat', 'w')
print>>fout, '#Index SIT  RA  Dec z kcorr mr v D LD Mr Mrk L'

for i in range(table.nrows):
    Index = table[0][i]
    SIT   = table[1][i]
    RA    = table[2][i]
    DEC   = table[3][i]
    z     = float(table[4][i])
    kcorr = float(table[5][i])
    mr    = float(table[6][i])

  
    v =c*z 
    #print v
    
    D = v/70.0
    
    LD = astCalc.dl(z)

    Mr = astCalc.absMag(mr,LD)
     
    Mrk = Mr-kcorr
    
    X = Ms-Mrk
    
    Y = 0.4*X
    
    L  = 10**Y
   

  
    
    print i+1,   Index, SIT, RA, DEC, z, kcorr, mr, v, D, LD, Mr, Mrk, L
  
    print>>fout, Index, SIT, RA, DEC, z, kcorr, mr, v, D, LD, Mr, Mrk, L
  
fout.close()


  
