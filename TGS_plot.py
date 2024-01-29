import os,asciidata
import numpy as np 
import math
import decimal
import matplotlib.pyplot as p
from astLib import astCoords
from astLib import astCalc
from scipy.stats.stats import pearsonr
from scipy import stats
from pylab import *

table = asciidata.open('TGS_trimatch.dat')

sigma_li=[] #sigma_1
rs_li =[]   #rs
rp_li=[]    #rp
Mvir_li=[]  #Mvir_rh
Mtrs_li=[]  # Mtrssun
Mtrsl_li=[] #Mtrsl
Mtrp_li=[]  #Mtrp
Mtrpl_li=[] #Mtrpl_3
rh_li=[]
rsF_li=[]
rsT_li=[]
l = (0,0.1,0.005)

f_out=open('Fsyst.dat','w')
print>>f_out,"#Index SIT RA DEC z  mr sigma  rs rp rh Mvirrh  Mtrs  Mtrsl Mtrp  Mtrpl rsF"  

fout= open('Tsyst.dat', 'w')
print>>fout, "#Index SIT RA DEC z  mr sigma  rs rp rh Mvirrh  Mtrs  Mtrsl Mtrp  Mtrpl rsT"

for i in range(table.nrows):
    Index_1 = int(table[0][i])
    SIT_1   = int(table[1][i])
    RA_1    = float(table[2][i])
    DEC_1   = float(table[3][i])
    z_1     = float(table[4][i])
    sigma_1 = float(table[5][i])
    rh      = float(table[6][i])
    Mvir_rh = float(table[7][i])
    R       = float(table[8][i])
    Mvir_rp = float(table[9][i])
    kcorr   = float(table[10][i])
    mr      = float(table[11][i])
    vmean   = float(table[12][i])
    sigmasq = float(table[13][i])
    rs12    = float(table[14][i])
    rs13    = float(table[15][i])
    rs23    = float(table[16][i])
    rs      = float(table[17][i])
    rpmean  = float(table[18][i])
    LS      = float(table[19][i])
    Reff    = float(table[20][i])
    Mtrssun = float(table[21][i])
    Mtrsl   = float(table[22][i])
    Mtrpsun = float(table[23][i])
    Mvir_rs = float(table[24][i])
    rp      = float(table[25][i])
    Mtrp    = float(table[26][i])
    Mtrpl_3 = float(table[27][i])
    rp12    = float(table[28][i])
    rp13    = float(table[29][i])
    rp23    = float(table[30][i])
    
    logrs= log(rs)
    print "logrs:", logrs
    
    sigma_li.append(sigma_1)
    rs_li.append(rs)
    rp_li.append(rp)
    Mvir_li.append(Mvir_rh)
    Mtrs_li.append(Mtrssun)
    Mtrsl_li.append(Mtrsl)
    Mtrp_li.append(Mtrp)
    Mtrpl_li.append(Mtrpl_3)
    rh_li.append(rh)
    
    sigma_array=np.asarray(sigma_li)
    #print "sigma_array:", sigma_array
    rs_array=np.asarray(rs_li)
    #print "rs_array:", rs_array
    rp_array=np.asarray(rp_li)
    Mvir_array=np.asarray(Mvir_li)
    Mtrp_array=np.asarray(Mtrp_li)
    Mtrpl_array=np.asarray(Mtrpl_li)
    Mtrs_array=np.asarray(Mtrs_li)
    Mtrsl_array=np.asarray(Mtrsl_li)
    rh_array=np.asarray(rh_li)
    
    
    if rs > 1.0:
	rsF=rs
	#print "rsF:", rsF
      
	rsF_li.append(rsF) 
	#print "Flist:", rsF_li
       
	FS=len(rsF_li)
	#print "FS:", FS
	
	print>>f_out, Index_1, SIT_1, RA_1, DEC_1, z_1,  mr, sigma_1,  rs, rp, rh, Mvir_rh,  Mtrssun,  Mtrsl, Mtrp,  Mtrpl_3, rsF  
    
   
    else:
	rsT=rs
	#print "rsT:", rsT
	
	rsT_li.append(rsT)
	#print "Tlist:", rsT_li
	
	TS=len(rsT_li)
	#print "TS:", TS
	
	print>>fout, Index_1, SIT_1, RA_1, DEC_1, z_1,  mr, sigma_1,  rs, rp, rh, Mvir_rh,  Mtrssun,  Mtrsl, Mtrp,  Mtrpl_3, rsT 



# to calculate correlation coeficient    

srs= pearsonr(sigma_array,rs_array)[0]
print "srs:",srs

srp=pearsonr(sigma_array,rp_array)[0]
print "srp:", srp

Mvrs=pearsonr(Mvir_array,rs_array)[0]
print "Mvrs:", Mvrs

Mvrp=pearsonr(Mvir_array,rp_array)[0]
print "Mvrp:", Mvrp

Mtlrp=pearsonr(Mtrpl_array,rp_array)[0]
print "Mtlrp:", Mtlrp

rsrp=pearsonr(rs_array,rp_array)[0]
print "rsrp:",rsrp

Mtlrs=pearsonr(Mtrsl_array,rs_array)[0]
print "Mtlrs:", Mtlrs

Mtrs=pearsonr(Mtrs_array,rs_array)[0]
print "Mtrs:", Mtrs

Mtrp=pearsonr(Mtrp_array,rp_array)[0]
print "Mtrp:", Mtrp
  
  
  
ms= mean(sigma_array)
print "meansigma:", ms

mds=median(sigma_array)
print "mediansigma:", mds
    
stds= np.std(sigma_li, ddof=1)
print "stdevsigma:", stds


mrh= mean(rh_array)
print "meanrh:", mrh

mdrh=median(rh_array)
print "medianrh:", mdrh
    
stdrh= np.std(rh_li, ddof=1)
print "stdevrh:", stdrh


mMvir= mean(Mvir_array)
print "meanMvir:", mMvir

mdMvir=median(Mvir_array)
print "medianMvir:", mMvir
    
stdMvir= np.std(Mvir_li, ddof=1)
print "stdevMvir:", stdMvir


mMtrp= mean(Mtrp_array)
print "meanMtrp:", mMtrp

mdMtrp=median(Mtrp_array)
print "medianMtrp:", mMtrp
    
stdMtrp= np.std(Mtrp_li, ddof=1)
print "stdevMtrp:", stdMtrp


mMtrs= mean(Mtrs_array)
print "meanMtrs:", mMtrs

mdMtrs=median(Mtrs_array)
print "medianMtrs:", mMtrs
    
stdMtrs= np.std(Mtrs_li, ddof=1)
print "stdevMtrs:", stdMtrs


mMtrpl= mean(Mtrpl_array)
print "meanMtrpl:", mMtrpl

mdMtrpl=median(Mtrpl_array)
print "medianMtrpl:", mMtrpl
    
stdMtrpl= np.std(Mtrpl_li, ddof=1)
print "stdevMtrpl:", stdMtrpl


mMtrsl= mean(Mtrsl_array)
print "meanMtrsl:", mMtrsl

mdMtrsl=median(Mtrsl_array)
print "medianMtrsl:", mMtrsl
    
stdMtrsl= np.std(Mtrsl_li, ddof=1)
print "stdevMtrsl:", stdMtrsl




## to plot these corr.coef.

p.figure(1)
p.plot([rs_li],[sigma_li],'ro')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('$r_s$ (Mpc)', fontsize=20)
p.ylabel ('$\sigma$ $(km s^{-1})$', fontsize=20)
#p.text(2.2, 125, r'(b)',fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot ([l],[l],'b-')
#p.show ()
savefig('rs_sigma.png')
os.system('convert rs_sigma.png rs_sigma.eps')
#os.system('cp rs_sigma.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(2)
p.plot([rp_li],[sigma_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('$\overline{r_p}$ (Mpc)', fontsize=20 )
p.ylabel ('$\sigma$ $(km s^{-1})$', fontsize=20)
#p.text(0.5, 120, r'(a)',fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot ([l],[l],'b-')
#p.subplot (1)
#p.show ()
savefig('rp_sigma.png')
os.system('convert rp_sigma.png rp_sigma.eps')
#copy the figure in the paper folder
os.system('cp rp_sigma.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

#p.figure(3)
#p.plot([rs_li],[rp_li],'ro')
#p.xlabel ('rs (Mpc)', fontsize=20)
#p.ylabel ('rp (Mpc)', fontsize=20)
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('rsrp.png')
#p.close()

#p.figure(4)
#p.plot([rs_li],[Mvir_li],'ro')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mvir/Msun')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('Mvirrs.png')
#p.close()

#p.figure(5)
#p.plot([rp_li],[Mvir_li],'r.')
#p.xlabel ('rp (Mpc)')
#p.ylabel ('Mvir/Msun')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('Mvirrp.png')
#p.close()

#p.figure(6)
#p.plot([rp_li],[Mtrp_li],'r.')
#p.xlabel ('rp (Mpc)')
#p.ylabel ('Mtrp/Msun')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('Mtrp.png')
#p.close()

#p.figure(7)
#p.plot([rs_li],[Mtrs_li],'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/Msun')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('Mtrs.png')
#p.close()

#p.figure(8)
#p.plot([rs_li],[Mtrsl_li],'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/L')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('MtrsL.png')
#p.close()

#p.figure(9)
#p.plot([rp_li],[Mtrpl_li],'r.')
#p.xlabel ('rp (Mpc)')
#p.ylabel ('Mtrp/L')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('MtrpL.png')
#p.close()

#p.figure(10)
#data= (Mvir_li)
#num_bins=25
#p.hist(data,num_bins)
#p.xlabel ('Mvir/Msun')
#p.ylabel ('Number of triplets')
##p.show ()
#savefig('Mvirrshisto.png')
#p.close()

p.figure(11)
data= (sigma_li)
num_bins=25
p.hist(data, num_bins, color='k', alpha=0.5)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$\sigma$  $(km s^{-1})$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(120, 25, r'(d)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('sigma_histo.png')
os.system('convert sigma_histo.png sigma_histo.eps')
#os.system('cp sigma_histo.eps ~/Amira_paper_I_v1/eps_figures/')
os.system('cp sigma_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(12)
data= (rh_li)
num_bins=25
p.hist(data, num_bins,color='k', alpha=0.5)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$r_h$ (Mpc)', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(0.5, 30, r'(c)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('rh_histo.png')
os.system('convert rh_histo.png rh_histo.eps')
#os.system('cp rh_histo.eps ~/Amira_paper_I_v1/eps_figures/')
os.system('cp rh_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(13)
data= (rs_li)
num_bins=25
p.hist(data, num_bins, color='k', alpha=0.5)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$r_s$ (Mpc)', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(2.5, 30, r'(b)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('rs_histo.png')
os.system('convert rs_histo.png rs_histo.eps')
#os.system('cp rs_histo.eps ~/Amira_paper_I_v1/eps_figures/')
os.system('cp rs_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

#p.figure(14)
#data= ( Mtrp_li)
#num_bins=25
#p.hist(data, num_bins)
#p.xlabel ('Mtrp/Msun')
#p.ylabel ('Number of triplets')
##p.show ()
#savefig('Mtrphisto.png')
#p.close()

#p.figure(15)
#p.plot([rh_li],[sigma_li],'r.')
#p.xlabel ('rh (Mpc)')
#p.ylabel ('Velocity Dispersion (km/s)')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('srh.png')
#p.close()

#p.figure(16)
#p.plot([rh_li],[Mtrpl_li],'r.')
#p.xlabel ('rh (Mpc)')
#p.ylabel ('Mtrp/L')
#p.plot ([l],[l],'b-')
##p.show ()
#savefig('Mtrh.png')
#p.close()

#p.figure(17)
#p.plot([rs_li],[Mtrp_li],'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/Msun')
#p.plot ([l],[l],'b-')
##sp.show ()
#savefig('Mtrpversusrs.png')
#p.close()


#x=rs_array
#y=sigma_array
#slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x, y)
#print "slope=", str(slope)
#print "r_value=", str(r_value)
#print "r_squared=", str(r_value**2)
#print "p_value=", str(p_value)
#predict_y= intercept+slope*x
##print predict_y
##pred_error=y-predict_y
##degrees_of_freedom= len(x)-2
##residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom   
#p.xlabel('rs(Mpc)')
#p.ylabel('Velocity dispersion (km/s)')
#p.plot(x,y,'o')
#p.plot(x,predict_y,'k-')
#p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
#savefig('equation_sigma_rs.jpg')
#p.close()



#p.figure(17)
#p.plot([rsT_li],[sigma_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Velocity Dispersion (km/s)')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[sigma_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Velocity Dispersion (km/s)')
#p.plot ([l],[l],'b-')
#savefig('srsTF.png')
#p.close()

#p.figure (18)
#p.plot([rsT_li],[Mvir_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mvir/Msun')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mvir_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mvir/Msun)')
#p.plot ([l],[l],'b-')
#savefig('MvirrsTF.png')
#p.close()

#p.figure (19)
#p.plot([rsT_li],[Mtrs_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/Msun')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mtrs_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/Msun)')
#p.plot ([l],[l],'b-')
#savefig('MtrsTF.png')
#p.close()

#p.figure (20)
#p.plot([rsT_li],[Mtrp_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/Msun')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mtrp_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/Msun)')
#p.figure (22)
#p.plot([rsT_li],[Mtrpl_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/L')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mtrpl_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/L)')
#p.plot ([l],[l],'b-')
#savefig('MtrplTF.png')
#p.close()p.plot ([l],[l],'b-')
#savefig('MtrpTF.png')
#p.close()


#p.figure (21)
#p.plot([rsT_li],[Mtrpl_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/L')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mtrpl_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrp/L)')
#p.plot ([l],[l],'b-')
#savefig('MtrplTF.png')
#p.close()


#p.figure (22)
#p.plot([rsT_li],[Mtrsl_li], 'r.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/L')
#p.plot ([l],[l],'b-')
#p.plot([rsF_li],[Mtrsl_li], 'b.')
#p.xlabel ('rs (Mpc)')
#p.ylabel ('Mtrs/L)')
#p.plot ([l],[l],'b-')
#savefig('MtrslTF.png')
#p.close()
