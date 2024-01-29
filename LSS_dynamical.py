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

table = asciidata.open('LSS_dynamical.dat')

sigma_li=[] #sigma_1
rs_li=[]
rp_li=[] #rp
Mvir_li=[]  #Mvir_rh
Mvirl_li=[]
d_li=[]
logd_li=[]
rh_li=[]
NLSS_li=[]
f_NLSS_li=[]
dNN_li=[]
KLSS_li=[]
QQLSS_li=[]
Qt_li=[]
Qratio_li=[]


for i in range(table.nrows):
    Index   = int(table[0][i])
    NLSS    = float(table[1][i])
    f_NLSS  = float(table[2][i])
    dNN     = float(table[3][i])
    KLSS    = float(table[4][i])
    QQLSS   = float(table[5][i])
    Qt      = float(table[6][i])
    Qratio  = float(table[7][i])
    sigma   = float(table[8][i])
    rh      = float(table[9][i])
    Mvir_rh = float(table[10][i])
    Mvirl   = float(table[11][i])
    rp      = float(table[12][i])
    v_mean  = float(table[13][i])
    d_mean  = float(table[14][i])
    rs      = float(table[15][i])
    
    
    #logd = log(d_mean)
    #print "logd:", logd
    
    #logMvir = log(Mvir_rh)
    
    #logMvirl = log(Mvirl)
   
   
   
    sigma_li.append(sigma)
    rs_li.append(rs)
    rp_li.append(rp)
    rh_li.append(rh)
    d_li.append(d_mean)
    Mvir_li.append(Mvir_rh)
    Mvirl_li.append(Mvirl)
    NLSS_li.append(NLSS)
    f_NLSS_li.append(f_NLSS)
    dNN_li.append(dNN)
    KLSS_li.append(KLSS)
    QQLSS_li.append(QQLSS)
    Qt_li.append(Qt)
    Qratio_li.append(Qratio)
    
    sigma_array=np.asarray(sigma_li)
    #print "sigma_array:", sigma_array
    d_array=np.asarray(d_li)
    dNN_array=np.asarray(dNN_li)
    rs_array=np.asarray(rs_li)
    rh_array=np.asarray(rh_li)
    rp_array=np.asarray(rp_li)
    Mvir_array=np.asarray(Mvir_li)
    Mvirl_array=np.asarray(Mvirl_li)
    NLSS_array=np.asarray(NLSS_li)
    f_NLSS_array=np.asarray(f_NLSS_li)
    KLSS_array=np.asarray(KLSS_li)
    QQLSS_array=np.asarray(QQLSS_li)
    Qt_array=np.asarray(Qt_li)
    Qratio_array=np.asarray(Qratio_li)
    
    
p.figure(1)
p.plot([dNN_li],[rh_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('dNN (Mpc)', fontsize=25)
p.ylabel ('$r_h$ (Mpc)', fontsize=25)
p.text(3, 0.5, r'(b)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('rh_dNN.png')
os.system('convert rh_dNN.png rh_dNN.eps')
os.system('cp rh_dNN.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(2)
p.plot([dNN_li],[sigma_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('dNN (Mpc)', fontsize=25)
p.ylabel ('$\sigma (km s^{-1})$', fontsize=25)
p.text(3, 120, r'(a)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('sigma_dNN.png')
os.system('convert sigma_dNN.png sigma_dNN.eps')
os.system('cp sigma_dNN.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(3)
p.plot([dNN_li],[Mvir_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('dNN (Mpc)', fontsize=25)
p.ylabel ('$M_{vir}/ M_{\odot} $', fontsize=25)
p.text(3, 1E13, r'(c)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvir_dNN.png')
os.system('convert Mvir_dNN.png Mvir_dNN.eps')
os.system('cp Mvir_dNN.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(4)
p.plot([dNN_li],[Mvirl_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('dNN (Mpc)', fontsize=25)
p.ylabel ('$M_{vir} / L_r$', fontsize=25)
p.text(3, 5000, r'(d)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvirl_dNN.png')
os.system('convert Mvirl_dNN.png Mvirl_dNN.eps')
os.system('cp Mvirl_dNN.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(5)
p.plot([QQLSS_li],[rh_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('QQLSS', fontsize=25)
p.ylabel ('$r_h$ (Mpc)', fontsize=25)
p.text(-3.5, 0.5, r'(b)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('rh_QQLSS.png')
os.system('convert rh_QQLSS.png rh_QQLSS.eps')
os.system('cp rh_QQLSS.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(6)
p.plot([QQLSS_li],[sigma_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('QQLSS', fontsize=25)
p.ylabel ('$\sigma (km s^{-1})$', fontsize=25)
p.text(-3.5, 120, r'(a)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('sigma_QQLSS.png')
os.system('convert sigma_QQLSS.png sigma_QQLSS.eps')
os.system('cp sigma_QQLSS.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(7)
p.plot([QQLSS_li],[Mvir_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('QQLSS', fontsize=25)
p.ylabel ('$M_{vir}/ M_{\odot} $', fontsize=25)
p.text(-3.5, 1E13, r'(c)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvir_QQLSS.png')
os.system('convert Mvir_QQLSS.png Mvir_QQLSS.eps')
os.system('cp Mvir_QQLSS.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(8)
p.plot([QQLSS_li],[Mvirl_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('QQLSS', fontsize=25)
p.ylabel ('$M_{vir}/ L_r$', fontsize=25)
p.text(-3.5, 5000, r'(d)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvirl_QQLSS.png')
os.system('convert Mvirl_QQLSS.png Mvirl_QQLSS.eps')
os.system('cp Mvirl_QQLSS.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

#p.figure(9)
#p.plot([Qt_li],[rh_li],'ro')
#p.rcParams.update({'xtick.labelsize': '15'})
#p.rcParams.update({'ytick.labelsize': '15'})
#p.xlabel ('Qt', fontsize=25)
#p.ylabel ('$r_h$ (Mpc)', fontsize=25)
##p.text(2.2, 0.5, r'(d)', fontsize=20) ## to insert a label inside the plot
#p.subplots_adjust(bottom=0.14, top=0.95)
##p.plot ([l],[l],'b-')
##p.show ()
#savefig('rh_Qt.png')
#os.system('convert rh_Qt.png rh_Qt.eps')
#os.system('cp rh_Qt.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
#p.close()

p.figure(9)
x=Qt_array
y=rh_array
slope, intercept, r_value, p_value, slope_std_error=stats.linregress(x,y)
print "slope=", str(slope)
print "r_value=", str(r_value)
print "r_squared=", str(r_value**2)
print "p_value=", str(p_value)
print "error=", str(slope_std_error)
predict_y= intercept+slope*x
#print predict_y
#pred_error=y-predict_y
#print "error=", pred_error
#degrees_of_freedom= len(x)-2
##residual_std_error=np.sqrt(np.sum(pred_error**2)/degrees_of_freedom
			   
p.rcParams.update({'axes.titlesize': '18'})
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('Qt', fontsize=25)
p.ylabel ('$r_h$ (Mpc)', fontsize=25)
p.text(1, 0.5, r'(b)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
p.plot(x,y,'ko')
p.xlim([-5,2])
p.plot(x,predict_y,'k-')
p.title('$y=%3.7sx+%3.7s$'%(slope, intercept))
savefig('rh_Qt.png')
os.system('convert rh_Qt.png rh_Qt.eps')
os.system('cp rh_Qt.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(10)
p.plot([Qt_li],[sigma_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('Qt', fontsize=25)
p.ylabel ('$\sigma (km s^{-1})$', fontsize=25)
p.text(2, 120, r'(a)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('sigma_Qt.png')
os.system('convert sigma_Qt.png sigma_Qt.eps')
os.system('cp sigma_Qt.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


p.figure(11)
p.plot([Qt_li],[Mvir_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('Qt', fontsize=25)
p.ylabel ('$M_{vir}/ M_{\odot} $', fontsize=25)
p.text(1.5, 1E13, r'(c)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvir_Qt.png')
os.system('convert Mvir_Qt.png Mvir_Qt.eps')
os.system('cp Mvir_Qt.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()

p.figure(12)
p.plot([Qt_li],[Mvirl_li],'ko')
p.rcParams.update({'xtick.labelsize': '15'})
p.rcParams.update({'ytick.labelsize': '15'})
p.xlabel ('Qt', fontsize=25)
p.ylabel ('$M_{vir}/ L_r$', fontsize=25)
p.text(2, 5000, r'(d)', fontsize=20) ## to insert a label inside the plot
p.subplots_adjust(bottom=0.14, top=0.95)
#p.plot ([l],[l],'b-')
#p.show ()
savefig('Mvirl_Qt.png')
os.system('convert Mvirl_Qt.png Mvirl_Qt.eps')
os.system('cp Mvirl_Qt.eps ~/2.Amira_newpaper/mnras_bold/eps_figures/')
p.close()


# to calculate correlation coeficient    

dNNsigma= pearsonr(dNN_array,sigma_array)[0]
print "dNNsigma:",dNNsigma

dNNrh= pearsonr(dNN_array,rh_array)[0]
print "dNNrh:",dNNrh

dNNrs= pearsonr(dNN_array,rs_array)[0]
print "dNNrs:",dNNrs

dNNrp=pearsonr(dNN_array,rp_array)[0]
print "dNNrp:", dNNrp

dNNMv=pearsonr(dNN_array, Mvir_array)[0]
print "dNNMv:", dNNMv

dNNMvl=pearsonr(dNN_array,Mvirl_array)[0]
print "dNNMvl:", dNNMvl

QQLSSsigma= pearsonr(QQLSS_array,sigma_array)[0]
print "QQLSSsigma:",QQLSSsigma

QQLSSrh= pearsonr(QQLSS_array,rh_array)[0]
print "QQLSSrh:",QQLSSrh

QQLSSrs= pearsonr(QQLSS_array,rs_array)[0]
print "QQLSSrs:",QQLSSrs

QQLSSrp=pearsonr(QQLSS_array,rp_array)[0]
print "QQLSSrp:", QQLSSrp

QQLSSMv=pearsonr(QQLSS_array, Mvir_array)[0]
print "QQLSSMv:", QQLSSMv

QQLSSMvl=pearsonr(QQLSS_array,Mvirl_array)[0]
print "QQLSSMvl:", QQLSSMvl

Qtsigma= pearsonr(Qt_array,sigma_array)[0]
print "Qtsigma:",Qtsigma

Qtrh= pearsonr(Qt_array,rh_array)[0]
print "Qtrh:",Qtrh

Qtrs= pearsonr(Qt_array,rs_array)[0]
print "Qtrs:",Qtrs

Qtrp=pearsonr(Qt_array,rp_array)[0]
print "Qtrp:", Qtrp

QtMv=pearsonr(Qt_array, Mvir_array)[0]
print "QtMv:", QtMv

QtMvl=pearsonr(Qt_array,Mvirl_array)[0]
print "QtMvl:", QtMvl