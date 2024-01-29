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

table = asciidata.open('dynamical_v_d.dat')

sigma_li=[] #sigma_1
rs_li=[]
rp_li=[] #rp
logMvir_li=[]  #logMvir_rh
logMvirl_li=[]
d_li=[]
logd_li=[]
rh_li=[]


for i in range(table.nrows):
    Index   = int(table[0][i])
    sigma   = int(table[1][i])
    rh      = float(table[2][i])
    Mvir_rh = float(table[3][i])
    Mvirl   = float(table[4][i])
    rp      = float(table[5][i])
    v_mean  = float(table[6][i])
    d_mean  = float(table[7][i])
    rs      = float(table[8][i])
    
    
    logd = log(d_mean)
    print "logd:", logd
    
    logMvir = log(Mvir_rh)
    
    logMvirl = log(Mvirl)
   
   
   
    sigma_li.append(sigma)
    rs_li.append(rs)
    rp_li.append(rp)
    rh_li.append(rh)
    d_li.append(d_mean)
    logd_li.append(logd)
    logMvir_li.append(logMvir)
    logMvirl_li.append(logMvirl)
   
    
    sigma_array=np.asarray(sigma_li)
    #print "sigma_array:", sigma_array
    d_array=np.asarray(d_li)
    logd_array=np.asarray(logd_li)
    rs_array=np.asarray(rs_li)
    rh_array=np.asarray(rh_li)
    rp_array=np.asarray(rp_li)
    logMvir_array=np.asarray(logMvir_li)
    logMvirl_array=np.asarray(logMvirl_li)
    
#figures with +ve correlations (1, 2, 3,4) (d vs sigma, Mvir, Mvir/L, rs)
p.figure(1)
x=d_array
y=sigma_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('$\overline{D}$ (Mpc)', fontsize=20)
p.ylabel('$\sigma (km s^{-1})$', fontsize=20)
p.text(0.5, 120, r'(a)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ko')
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('sigma_d.png')
os.system('convert sigma_d.png sigma_d.eps')
os.system('cp sigma_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(2)
x=logd_array
y=logMvir_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('log($\overline{D}$)', fontsize=20)
p.ylabel('log($M_{vir}$)', fontsize=20)
p.text(-2.0, 30, r'(b)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ko')
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('Mvir_d.png')
os.system('convert Mvir_d.png Mvir_d.eps')
os.system('cp Mvir_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(3)
x=logd_array
y=logMvirl_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('log($\overline{D}$)', fontsize=20)
p.ylabel('log($M_{vir}/L_r$)', fontsize=20)
p.text(-2.0, 8, r'(c)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ko')
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('Mvirl_d.png')
os.system('convert Mvirl_d.png Mvirl_d.eps')
os.system('cp Mvirl_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(4)
x=d_array
y=rs_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#degrees_of_freedom= len(x)-2
#residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel('$\overline{D}$ (Mpc)', fontsize=20)
p.ylabel('$\overline{r_s}$ (Mpc)', fontsize=20)
p.text(0.5, 2.5, r'(f)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ko')
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('rs_d.png')
os.system('convert rs_d.png rs_d.eps')
os.system('cp rs_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

# figures with scattered relations (5, 6) (d vs rh, rp)
p.figure(5)
p.plot([d_li],[rh_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('$\overline{D}$ (Mpc)', fontsize=25)
p.ylabel ('$r_h$ (Mpc)', fontsize=25)
p.text(2.2, 0.5, r'(d)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('rh_d.png')
os.system('convert rh_d.png rh_d.eps')
os.system('cp rh_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(6)
p.plot([d_li],[rp_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('$\overline{D}$ (Mpc)', fontsize=25)
p.ylabel ('$\overline{r_p}$ (Mpc)', fontsize=25)
p.text(2.2, 0.5, r'(e)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('rp_d.png')
os.system('convert rp_d.png rp_d.eps')
os.system('cp rp_d.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


# to calculate correlation coeficient    

dsigma= pearsonr(d_array,sigma_array)[0]
print "dsigma:",dsigma

drh= pearsonr(d_array,rh_array)[0]
print "drh:",drh

drs= pearsonr(d_array,rs_array)[0]
print "drs:",drs

drp=pearsonr(d_array,rp_array)[0]
print "drp:", drp

dMv=pearsonr(logd_array, logMvir_array)[0]
print "dMv:", dMv

dMvl=pearsonr(logd_array,logMvirl_array)[0]
print "dMvl:", dMvl