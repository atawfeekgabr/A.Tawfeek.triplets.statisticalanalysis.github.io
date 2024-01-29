#to read table
import os,asciidata
import numpy as np
from astLib import astCoords
from astLib import astCalc

table = asciidata.open('SIT_input5_new.dat')

fout = open('TGS.dat', 'w')
print>>fout, '#Index SIT RA DEC z kcorr mr'

#creat an array in a new table from 1 to 315   
#ind_array = np.arange(1, 316)

#print ind_array
#print>>fout,"#ind_array"

for i in range(table.nrows):
    Index = table[0][i]
    SIT   = table[1][i]
    RA    = table[2][i]
    DEC   = table[3][i]
    z     = float(table[4][i])
    kcorr = float(table[5][i])
    mr    = float(table[6][i])
    
    
    if (SIT==1):
      
      print i+1,   Index, SIT, RA, DEC, z, kcorr, mr
      
      print>>fout, Index, SIT, RA, DEC, z, kcorr, mr

#creat an array in a new table from 1 to 315   
#ind_array = np.arange(1, 316)

#print ind_array

#for item in len(ind_array)):

#print item

#print>>fout, item


fout.close()
