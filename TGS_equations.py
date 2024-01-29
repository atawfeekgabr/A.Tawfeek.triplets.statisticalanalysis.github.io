import os,asciidata
import numpy as np 
import math
import decimal
import matplotlib.pyplot as p
from matplotlib.ticker import MaxNLocator
from astLib import astCoords
from astLib import astCalc
from scipy.stats.stats import pearsonr
from scipy import stats
from pylab import *

table = asciidata.open('logMvir_rh_rs.dat')

sigma_li=[] #sigma_1
logrs_li =[]   #logrs
rs_li=[]
rp_li=[]    #rp
logMvir_li=[]  #logMvir_rh
logMvirl_li=[]
logrp_li=[]

for i in range(table.nrows):
    Index_1 = int(table[0][i])
    SIT_1   = int(table[1][i])
    RA_1    = float(table[2][i])
    DEC_1   = float(table[3][i])
    z_1     = float(table[4][i])
    sigma_1 = float(table[5][i])
    rh      = float(table[6][i])
    Mvir_rh = float(table[7][i])
    logMvir = float(table[8][i])
    Mvirl   = float(table[9][i])
    logMvirl= float(table[10][i])
    rp      = float(table[11][i])
    rs      = float(table[12][i])
    logrs   = float(table[13][i])
    
    logrp=log(rp)
    print "logrp:", logrp
   
   
   
    sigma_li.append(sigma_1)
    rs_li.append(rs)
    logrs_li.append(logrs)
    rp_li.append(rp)
    logMvir_li.append(logMvir)
    logMvirl_li.append(logMvirl)
    logrp_li.append(logrp)
    
    sigma_array=np.asarray(sigma_li)
    #print "sigma_array:", sigma_array
    logrs_array=np.asarray(logrs_li)
    #print "rs_array:", rs_array
    rp_array=np.asarray(rp_li)
    rs_array=np.asarray(rs_li)
    logMvir_array=np.asarray(logMvir_li)
    logMvirl_array=np.asarray(logMvirl_li)
    
    
p.figure(1)
x=logrs_array
y=logMvir_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
print "error=", str(slope_std_error)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('log($\overline{r_s}$)', fontsize=20)
p.ylabel('log($M_{vir}$)', fontsize=20)
p.text(-0.5, 13.0, r'(b)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ro')
p.xlim([-0.8,0.48])
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('eqn_Mvir_rs.png')
os.system('convert eqn_Mvir_rs.png eqn_Mvir_rs.eps')
os.system('cp eqn_Mvir_rs.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()




p.figure(2)
x=logrs_array
y=logMvirl_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
print "error=", str(slope_std_error)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
#p.rcParams.update({'legend.fontsize': 'x-large'})
#p.rcParams.update({'axes.labelsize': 'x-large'})
p.xlabel('log($\overline{r_s}$)', fontsize=20)
p.ylabel('log($M_{vir}/L_r$)', fontsize=20)
p.text(-0.5, 3.6, r'(b)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
p.subplots_adjust(hspace=.5)
p.subplots_adjust(wspace=2.8)
p.plot(x,y,'ro')
p.xlim([-0.8,0.48])
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('eqn_ML_rs.png')
os.system('convert eqn_ML_rs.png eqn_ML_rs.eps')
os.system('cp eqn_ML_rs.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()



p.figure(3)
x=rs_array
y=sigma_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x, y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
print "error=", str(slope_std_error)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('$\overline{r_s}$ (Mpc)', fontsize=20)
p.ylabel('$\sigma$  $(km s^{-1})$', fontsize=20)
p.text(0.5, 120, r'(b)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
p.subplots_adjust(hspace=.5)
p.plot(x,y,'ro')
p.xlim([0.2,2.8])
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx%3.7s$'%(slope, intercept))
savefig('eqn_rs_sigma.png')
os.system('convert eqn_rs_sigma.png eqn_rs_sigma.eps')
os.system('cp eqn_rs_sigma.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()
#p.show()

p.figure(4)
p.plot([logrs_li],[logMvir_li],'ro')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('log($\overline{r_s}$)', fontsize=25)
p.ylabel ('log($M_{vir}$)', fontsize=25)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvir_rs.png')
os.system('convert Mvir_rs.png Mvir_rs.eps')
os.system('cp Mvir_rs.eps ~/Amira_paper_I_v1/eps_figures/')
p.close()

p.figure(5)
p.plot([logrs_li],[logMvirl_li],'ro')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('log($\overline{r_s}$)', fontsize=40)
p.ylabel ('log($M_{vir}/L_r$)', fontsize=40)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('ML_rs.png')
os.system('convert ML_rs.png ML_rs.eps')
os.system('cp ML_rs.eps ~/Amira_paper_I_v1/eps_figures/')
p.close()

p.figure(6)
p.plot([logrp_li],[logMvir_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.rcParams.update({'legend.fontsize': 'x-large'})
p.rcParams.update({'axes.labelsize': 'x-large'})
p.xlabel ('log($\overline{r_p}$)', fontsize=20)
p.ylabel ('log($M_{vir}$)', fontsize=20)
#p.text(-3.0, 13, r'(a)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvir_rp.png')
os.system('convert Mvir_rp.png Mvir_rp.eps')
os.system('cp Mvir_rp.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(7)
p.plot([logrp_li],[logMvirl_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('log($\overline{r_p}$)', fontsize=20)
p.ylabel ('log($M_{vir}/L_r$)', fontsize=20)
#p.text(-2.8, 3.7, r'(a)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('ML_rp.png')
os.system('convert ML_rp.png ML_rp.eps')
os.system('cp ML_rp.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()