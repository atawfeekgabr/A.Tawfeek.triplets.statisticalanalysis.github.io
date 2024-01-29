#to read table
from astLib import astCoords
from astLib import astCalc
from numpy import *
from pylab import *

file = "SIT_input2.dat"
d =loadtxt(file)

c=2.99792*10**5
#PIG=9pi/2G (constant for virial mass equation)
PIG=2.118*10**11
# solar absolute mag in r-bad
Ms=4.68
#print c

fout = open('SIT_newoutput.dat', 'w')
print>>fout, '#Index SIT  RA  Dec z kcorr mr v D LD Mr Mrk  L '

for i in range(d.nrows):
    Index = int(d[0][i])
    SIT   = int(d[1][i])
    RA    = float(d[2][i])
    DEC   = float(d[3][i])
    z     = float(d[4][i])
    kcorr = float(d[5][i])
    mr    = float(d[6][i])

  
    v =c*z 
    print v
    D = v/70.0
    
    LD = astCalc.dl(z)

    Mr = astCalc.absMag(mr,LD)
     
    Mrk = Mr-kcorr
    
    X = Ms-Mrk
    
    Y = 0.4*X
    
    L  = 10**Y
   

  
    
    print i+1, Index, SIT, RA, DEC, z, kcorr, mr, v, D, LD, Mr, Mrk, L
  
    print>>fout, Index, SIT, RA, DEC, z,kcorr, mr, v, D, LD, Mr, Mrk, L
  
fout.close()

#os.system('cp SIT_output3.dat SIT_new2.dat')
  
