#to filter the galaxy table  

import os, asciidata
import matplotlib.pyplot as p
from astLib import astCoords
from pylab import *
from scipy.stats.stats import pearsonr

# open the sdsss data table.
table1 = asciidata.open('newoutput_AmiraGabr.dat')


#create a new file for the output
fout = open('CasJobfz.dat', 'w')
print>>fout, "#id_CLG  ra_CLG   dec_CLG   objID	ra	dec	g	err_g	r	err_r	  i	err_i	  z	err_z	zp	err_zp	 kcorrR lumDist Mr  zs	  err_zs  velDisp class   sep_arcmin   sep_arcsec  "

objID_li = []
zp_li = []
zs_li = []	
l = (0,0.1,0.005) 

print l


# to loop over the table1
for k in range (table1.nrows):
    id_CLG        = table1[0][k]
    ra_CLG        = table1[1][k]
    dec_CLG       = table1[2][k]
    objID         = table1[3][k]
    ra            = table1[4][k]
    dec           = table1[5][k]
    g             = float(table1[6][k])
    err_g         = float(table1[7][k])
    r             = float(table1[8][k])
    err_r         = float(table1[9][k])
    i             = float(table1[10][k])
    err_i         = float(table1[11][k])
    z             = float(table1[12][k])
    err_z         = float(table1[13][k])    
    zp            = float(table1[14][k])
    err_zp        = float(table1[15][k])
    kcorrR        = float(table1[16][k])
    lumDist       = float(table1[17][k])
    Mr            = float(table1[18][k])
    zs            = float(table1[19][k])
    err_zs        = float(table1[20][k])
    velDisp       = float(table1[21][k])
    gclass        = str(table1[22][k])
                   
        
       
#select galaxies with photometric redshift
    if zp > 0.0 and zs > 0.0 :
      zp_li.append(zp)
      zs_li.append(zs)
      
zp_array=np.asarray(zp_li)
zs_array=np.asarray(zs_li)
      
zpzs= pearsonr(zp_array,zs_array)[0]
print "zpzs:", zpzs
      
      
      #print len(zp_li), len(zs_li)
      
      
p.plot([zp_li],[zs_li],'r.')
p.ylabel ('zs')
p.xlabel ('zp')
p.plot ([l],[l],'b-')
p.show ()
# select galaxies that are classified as galaxies from their spectra
        #if gclass in ['GALAXY', '-9999'] :
                    
#to avoid repeating sources    
            #if long(objID) not in objID_li:
                        
                #objID_li.append(long(objID))
			    
 #calculate the angular separation in degrees, arcmin and arcsec.
                #sep_deg    =  astCoords.calcAngSepDeg(ra_CLG, dec_CLG, ra, dec)
                #sep_arcmin =  float(sep_deg*60.0) 
                #sep_arcsec =  float(sep_deg*3600)               
                                            
                #print  k+1,   long(id_CLG),  ra_CLG,  dec_CLG,  long(objID),  ra,  dec,  g,  err_g,  r,  err_r,  i,  err_i, z, err_z,  zp,  err_zp, kcorrR, lumDist, Mr, zs, err_zs, velDisp, gclass
                #print>>fout, long(id_CLG),  ra_CLG,  dec_CLG,  long(objID),  ra,  dec,  g,  err_g,  r,  err_r,  i,  err_i, z, err_z,  zp,  err_zp,  kcorrR, lumDist, Mr, zs, err_zs, velDisp, gclass
 # close the output file 
fout.close()

 #number of the output galaxies
#print "len(objID_li):", len(objID_li)

