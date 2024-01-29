import os,asciidata
import numpy as np
import math
import decimal
import matplotlib.pyplot as p
from astLib import astCoords
from astLib import astCalc
from scipy.stats.stats import pearsonr
#import statistics
from pylab import *

c=2.99792*10**5
PIG=2.118*10**11
Ms=4.68
MPC_km= 3.0857*10**19 #from MPC to meters

M1_li=[]
M2_li=[]
M3_li=[]

L1_li=[] 
L2_li=[]
L3_li=[]


M12_li=[]
M13_li=[]
M1M23_li=[]

L12_li=[]
L13_li=[]
L1L23_li=[]

table= asciidata.open('SIT_123.dat')
fout = open('SIT_M_Lratios.dat', 'w')
f_out= open ('SIT_ML1_ML23.dat', 'w')
print>>fout, "#Index M1/M2 M1/M3 M12_M13 L1/L2 L1/L3 L12_L13"
print>>f_out, "#Index M1/M23 L1/L23"

for i in range(table.nrows):
    Index = int(table[0][i])
    RA1   = float(table[1][i])
    DEC1  = float(table[2][i])
    RA2   = float(table[3][i])
    DEC2  = float(table[4][i])
    RA3   = float(table[5][i])
    DEC3  = float(table[6][i])
    z1    = float(table[7][i])
    z2    = float(table[8][i])
    z3    = float(table[9][i])
    v1    = float(table[10][i])
    v2    = float(table[11][i])
    v3    = float(table[12][i])
    d1    = float(table[13][i])
    d2    = float(table[14][i])
    d3    = float(table[15][i])
    Mrk1  = float(table[16][i])
    Mrk2  = float(table[17][i])
    Mrk3  = float(table[18][i])
    L1    = float(table[19][i])
    L2    = float(table[20][i])
    L3    = float(table[21][i])
    
    
    
    M1_M2 = Mrk1/Mrk2
    M1_M3 = Mrk1/Mrk3
    M12_M13 = abs( (M1_M2) - (M1_M3))
    print "M12-M13:", M12_M13
    
    M1_23 = Mrk1/(Mrk2+Mrk3)
    
    L1_L2= L1/L2
    L1_L3= L1/L3
    L12_L13= abs((L1_L2) - (L1_L3))
    print "L12-L13:", L12_L13
    
    L1_23= L1/(L2+L3)
    
    M12_li.append(M1_M2)
    M13_li.append(M1_M3)
    M1M23_li.append(M1_23)
    
    L12_li.append(L1_L2)
    L13_li.append(L1_L3)
    L1L23_li.append(L1_23)
    
    M1_li.append(Mrk1)
    M2_li.append(Mrk2)
    M3_li.append(Mrk3)
    
    L1_li.append(L1)
    L2_li.append(L2)
    L3_li.append(L3)
    
    mn_M1=np. mean(M1_li)
    print "M1mean:", mn_M1
    mn_M2=np. mean(M2_li)
    print "M2mean:", mn_M2
    mn_M3=np. mean(M3_li)
    print "M3mean:", mn_M3
    mn_M12=np. mean(M12_li)
    print "M12mean:", mn_M12
    mn_M13=np. mean(M13_li)
    print "M13mean:", mn_M13
    
    mn_L1=np. mean(L1_li)
    print "L1mean:", mn_L1
    mn_L2=np. mean(L2_li)
    print "L2mean:", mn_L2
    mn_L3=np. mean(L3_li)
    print "L3mean:", mn_L3
    mn_L12=np. mean(L12_li)
    print "L12mean:", mn_L12
    mn_L13=np. mean(L13_li)
    print "L13mean:", mn_L13
    
    
    md_M1=np. median(M1_li)
    print "M1median:", md_M1
    md_M2=np. median(M2_li)
    print "M2median:", md_M2
    md_M3=np. median(M3_li)
    print "M3median:", md_M3
    md_M12=np. median(M12_li)
    print "M12median:", md_M12
    md_M13=np. median(M13_li)
    print "M31median:", md_M13
    
    md_L1=np. median(L1_li)
    print "L1median:", md_L1
    md_L2=np. median(L2_li)
    print "L2median:", md_L2
    md_L3=np. median(L3_li)
    print "L3median:", md_L3
    md_L12=np. median(L12_li)
    print "L12median:", md_L12
    md_L13=np. median(L13_li)
    print "L13median:", md_L13
    
    std_M1=np.std(M1_li)
    print "M1std:", std_M1
    std_M2=np. std(M2_li)
    print "M2std:", std_M2
    std_M3=np. std(M3_li)
    print "M3std:", std_M3
    std_M12=np. std(M12_li)
    print "M12std:", std_M12
    std_M13=np. std(M13_li)
    print "M13std:", std_M13
    
    std_L1=np. std(L1_li)
    print "L1std:", std_L1
    std_L2=np. std(L2_li)
    print "L2std:", std_L2
    std_L3=np.std(L3_li)
    print "L3std:", std_L3
    std_L12=np. std(L12_li)
    print "L12std:", std_L12
    std_L13=np.std(L13_li)
    print "L13std:", std_L13
    
    
    print>>fout, Index, M1_M2, M1_M3, M12_M13, L1_L2, L1_L3, L12_L13
    print>>f_out, Index, M1_23, L1_23
    
    
    
p.figure(1)
data= ( M12_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$M_{rG1}/M_{rG2}$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(1.30, 35, r'(a)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('M12_histo.png')
os.system('convert M12_histo.png M12_histo.eps')
os.system('cp M12_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(2)
data= ( M13_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$M_{rG1}/M_{rG3}$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(1.45, 50, r'(b)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('M13_histo.png')
os.system('convert M13_histo.png M13_histo.eps')
os.system('cp M13_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(3)
data= ( M1M23_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$M_{rG1}/(M_{rG2}+M_{rG3})$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(0.66, 35, r'(c)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('M1M23_histo.png')
os.system('convert M1M23_histo.png M1M23_histo.eps')
os.system('cp M1M23_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(4)
data= ( L12_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$L_{rG1}/L_{rG2}$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(80, 170, r'(d)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('L12_histo.png')
os.system('convert L12_histo.png L12_histo.eps')
os.system('cp L12_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(5)
data= ( L13_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$L_{rG1}/L_{rG3}$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(600, 255, r'(e)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('L13_histo.png')
os.system('convert L13_histo.png L13_histo.eps')
os.system('cp L13_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()

p.figure(6)
data= ( L1L23_li)
num_bins=25
p.hist(data, num_bins)
p.rcParams.update({'xtick.labelsize': '13'})
p.rcParams.update({'ytick.labelsize': '13'})
p.xlabel ('$L_{rG1}/(L_{rG2}+L_{rG3})$', fontsize=25)
p.ylabel ('Number of triplets', fontsize=25)
p.text(40, 140, r'(f)', fontsize=20)
p.subplots_adjust(bottom=0.14, top=0.95)
#p.show ()
savefig('L1L23_histo.png')
os.system('convert L1L23_histo.png L1L23_histo.eps')
os.system('cp L1L23_histo.eps ~/Amira_paper_I_v1/MNRAS_ApSS_templates/mnras/eps_figures/')
p.close()