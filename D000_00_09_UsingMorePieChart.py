'''
Created on Apr 16, 2019
Beginning DDA Modeling
Kei215-4~7
@author: keiin
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *
from af002_energy_calculation import *
from N2fixBudget import *
from E_BiosynthesisFromNH4 import *

def QcFromV(V): #Compute Qc from V #Menden-Deuer 2000
    return 0.216*V**0.939/12 #(pmol C cell-1)

def QcFromVdiatom(V):
    return 0.288*V**0.811/12 #(pmol C cell-1)

def QnFromV(V): #Compute Qn from V #Menden-Deuer 2000
    return 0.118*V**0.849/14 #(pmol N cell-1)

CN = 6.6

TrichoData = genfromtxt('..\\Data\\TrichoMuN2fix.csv',delimiter=',').T
CrocoData = genfromtxt('..\\Data\\CrocoMuN2fix.csv',delimiter=',').T
NostocData = genfromtxt('..\\Data\\NostocMuN2fix.csv',delimiter=',').T

Nv = 4 #Number of the vegetative cells per heterocysts
Nh = 2 #Number of the heterocysts per host diatom

Vv0 = 18.845 #(um3) Cell Volume of the vegetative cells from "03 Cell volume.xlsx"
Vh0 = 61.0125 #(um3) Cell volume of the heterocysts "03 Cell volume.xlsx"
Vd0 = 3493.529 #(um3) Cell volume of the host diatom "03 Cell volume.xlsx"

Vv = Vv0*Nv*Nh #(um3) Cell Volume of the vegetative cells from "03 Cell volume.xlsx"
Vh = Vh0*Nh #(um3) Cell volume of the heterocysts "03 Cell volume.xlsx"
Vd = Vd0 #(um3) Cell volume of the host diatom "03 Cell volume.xlsx"

QcV0 = QcFromV(Vv0) #(pmol C cell-1)
QcH0 = QcFromV(Vh0) #(pmol C cell-1)
QcD0 = QcFromVdiatom(Vd0) #(pmol C cell-1)

QcV = QcV0*Nv*Nh
QcH = QcH0*Nh
QcD = QcD0

QcC = QcV + QcH #C as "Chain"
Qc = QcV + QcH + QcD

print(QcH/QcC)

QnV = QcV/CN #(pmol N cell-1)
QnH = QcH/CN #(pmol N cell-1)
QnD = QcD/CN #(pmol N cell-1)

QnC = QnV + QnH
Qn = QnV + QnH + QnD

print(QcV/QnV,QcH/QnH,QcD/QnD)

#Maybe I can make it growth rate array and use it for the x-axis.
MuMax = 0.77 #(d-1) Maximum growth rate; Villareal 1989 (from Foster 2010)
MuStep = 0.01
Mu = arange(0.0,MuMax+MuStep,MuStep) #(d-1) Growth rate

Fn2fixN = Mu*(QnV + QnH + QnD) #(pmol N cell-1 d-1) 

E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate

YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron

print(YcnN2fix,YresN2fix)

Fn2fixC = Fn2fixN*YcnN2fix + Fn2fixN*YresN2fix #(pmol C cell-1 d-1) C consumption for N2 fixation

Fpho = Mu*(QcV + QcH + QcD)*(1+E) + Fn2fixC #(pmol C cell-1 d-1) Photosynthesis rate

FphoV = Fpho*QnV/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by vegetative cells
FphoD = Fpho*QnD/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by diatom cells

Xlabel = '$\mathit{\mu}$ (d$^{-1}$)'

PanelPlot = 0
if PanelPlot == 1:
    figure(1,figsize=(25,12))
    rcParams.update({'font.size': 15})
    
def sup(n):
    if PanelPlot == 1:
        subplot(2,4,n)
    else:
        figure(n)
    xlabel(Xlabel)

def sf(name):
    if PanelPlot != 1:
        Savefig3('02\\04 DDA',name,300)

def ti(name):
    title(name,y=1.02)

sup(1)
stackplot(Mu,FphoV,FphoD,labels=['Vegetative','Diatom'])
ylabel('pmol C cell$^{-1}$ d$^{-1}$')
ti('C supply (Photosynthesis)')
legend(loc = 2)
sf('PhotoSynthesis')

sup(2)
plot(Mu,Fn2fixN/QcC,label='Model, DDA',color="#1F77B4")
plot(TrichoData[0],TrichoData[1],'o',label='Data, Tricho',color="#FF7F0E")
plot(CrocoData[0],CrocoData[1],'o',label='Data, Croco',color="#2CA02C")
plot(NostocData[0],NostocData[1],'o',label='Data, Nostoc',color="#9467BD")
plot(Mu,Fn2fixN/QcC*QnC/Qn,'--',label='Model, DDA*',color="k")#D62728")

ylim(top=0.75)
ylabel('mol N mol C$^{-1}$ d$^{-1}$')
ti('N$_2$ fixation')
legend(loc = 2,fontsize=17)
sf('N2fix')

sup(3)
stackplot(Mu,Mu*QnV/QcC,Mu*QnH/QcC,Mu*QnD/QcC,labels=['Vegetative','Heterocyst','Diatom'],colors=['#1F77B4','#2CA02C','#FF7F0E'])
plot(Mu,Fn2fixN/QcC,label='Model, DDA',color="#1F77B4")
plot(Mu,Fn2fixN/QcC*QnC/Qn,'--',label='Model, DDA*',color="k")
ti('Fate of N')
ylabel('mol N mol C$^{-1}$ d$^{-1}$')
legend(loc = 2)
sf('Fate N')

sup(8)
plot(Mu,FphoD,'--',color='#FF7F0E',label='Supp. Diatom')
plot(Mu,Mu*QcD*(1+E),color='#FF7F0E',label='Cons. Diatom')
plot(Mu,Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC,color='#1F77B4',label='Cons. Chain')
plot(Mu,FphoV,'--',color='#1F77B4',label='Supp. Chain')
legend(loc = 2)
ylabel('pmol C cell$^{-1}$ d$^{-1}$')
ti('C Supply and Consumption')
sf('C_Supply_Consumption')

figure(9)
plot(Mu,Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC,':',color="#2CA02C",label="Consumption")
plot(Mu,Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC*(QnV+QnH)/Qn,color="#2CA02C",label="Consumption*")
plot(Mu,FphoV,'--',color="#2CA02C",label="Supply")
xlabel(Xlabel)
ti("C budget for Chain")
ylabel('pmol C cell$^{-1}$ d$^{-1}$')
legend(loc = 2)
sf("C_budget_chain")

#+++++++++++++++++++++++++++++++++++++
# Plotting Fate of C in DDA
#+++++++++++++++++++++++++++++++++++++
a = Mu*(QcD+QcV+QcH)/Fpho*100
b = Mu*(QcD+QcV+QcH)*E/Fpho*100
c = Fn2fixC/Fpho*100
values = [nanmean(a),nanmean(b),nanmean(c)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
if round(c[-1],10) != round(nanmean(c),10):
    print("no ,good")
#-----------------------------

figure(12)
pie(values,labels=['Biomass','Syn. Resp.','N$_{2}$ fixation'],colors=['#F4B084','#8EA9DB','#A9D08E'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Fate of C in DDA")
sf("Pie_Fate_of_C_in_DDA")

#+++++++++++++++++++++++++++++++++++++
# Plotting Fate of C fixed by diatom
#+++++++++++++++++++++++++++++++++++++
Ctrans = ((Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC)-FphoV)/FphoD*100
a = Ctrans
b = 100 - Ctrans
values = [nanmean(a),nanmean(b)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
#-----------------------------

figure(13)
pie(values,labels=['Chain','Diatom'],colors=['#00B050','#FFC000'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Fate of C fixed by diatom")
sf("Pie_Fate_C_from_Diatom")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# Plotting Fate of N fixed by the Chain
#++++++++++++++++++++++++++++++++++++++++

a = Mu*QnV/QcC/Fn2fixN*100
b = Mu*QnH/QcC/Fn2fixN*100
c = Mu*QnD/QcC/Fn2fixN*100

values = [nanmean(a),nanmean(b),nanmean(c)]

#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
def Print(x):
    print(round(x[-1],10),round(nanmean(x),10))
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)
if round(c[-1],10) != round(nanmean(c),10):
    print("no good")
    Print(c)
#-----------------------------

figure(14)
pie(values,labels=['Veg','Het','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Fate of N")
sf("Pie_Fate_N")
N = nanmean(c)/(nanmean(a)+nanmean(b))
print('N',N)
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C source in the chain
#++++++++++++++++++++++++++++++++++++++++
Cost = Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC
a = FphoV/Cost*100
b = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC-FphoV)/Cost*100

values = [nanmean(a),nanmean(b)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)

#-----------------------------

figure(18)
pie(values,labels=['Veg','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C source in the Chain")
sf("Pie_C_source_Chain")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C source in DDA
#++++++++++++++++++++++++++++++++++++++++

a = FphoV/Fpho*100
b = FphoD/Fpho*100

values = [nanmean(a),nanmean(b)]
print(values)
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)

#-----------------------------

figure(15)
pie(values,labels=['Veg','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C source in DDA")
sf("Pie_C_source_DDA")
#++++++++++++++++++++++++++++++++++++++

if PanelPlot == 1:
    Savefig3('02\\04 DDA','LargePlot',300)
    
#++++++++++++++++++++++++++++++++++++++++
# C cost in DDA
#++++++++++++++++++++++++++++++++++++++++
tot = Mu*QcD*(1+E) + Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC
a = Mu*QcD*(1+E)/tot*100
b = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC)/tot*100

values = [nanmean(a),nanmean(b)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)

#-----------------------------

figure(16)
pie(values,labels=['Diatom','Chain'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C cost in DDA")
sf("Pie_C_cost_DDA")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C source for the chain
#++++++++++++++++++++++++++++++++++++++++
tot = Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC*(QnV+QnH)/Qn
a = FphoV/tot*100
b = (tot-FphoV)/tot*100

values = [nanmean(a),nanmean(b)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)

#-----------------------------
figure(17)
pie(values,labels=['Chain','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C source in Chain")
sf("Pie_C_budget_Chain")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# Volume
#++++++++++++++++++++++++++++++++++++++++
V = Vv + Vh + Vd
a = Vv/V*100
b = Vh/V*100
c = Vd/V*100

values = [a,b,c]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
#-----------------------------
figure(19)
pie(values,labels=['Veg','Het','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Volume")
sf("Pie_Volume")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C in Chain
#++++++++++++++++++++++++++++++++++++++++
QcC = QcV + QcH
a = QcV/QcC*100
b = QcH/QcC*100

values = [a,b]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
#-----------------------------
figure(21)
pie(values,labels=['Veg','Het'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Q$_{C}$")
sf("Pie_Qc")
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# N in DDA
#++++++++++++++++++++++++++++++++++++++++
Qn = QnV + QnH + QnD
a = QnV/Qn*100
b = QnH/Qn*100
c = QnD/Qn*100

values = [a,b,c]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
#-----------------------------
figure(22)
pie(values,labels=['Veg','Het','Diatom'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Q$_{N}$")
sf("Pie_Qn")
print('Qn',c/(a+b))
#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C consumption Chain
#++++++++++++++++++++++++++++++++++++++++
tot = Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC

a = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E)/tot*100
b = Fn2fixC/tot*100

values = [nanmean(a),nanmean(b)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)
#-----------------------------

figure(23)
pie(values,labels=['Veg','Het'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C consumption in the Chain")
sf("Pie_C_consumpt_Chain")

#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# C consumption in DDA
#++++++++++++++++++++++++++++++++++++++++
tot = Mu*QcD*(1+E)+Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC
a = Mu*QcD*(1+E)/tot*100
b = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E)/tot*100
c = Fn2fixC/tot*100

values = [nanmean(a),nanmean(b),nanmean(c)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)
if round(c[-1],10) != round(nanmean(c),10):
    print("no good")
    Print(c)
#-----------------------------

figure(24)
pie(values,labels=['Dia','Veg','Het'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("C consumption in DDA")
sf("Pie_C_consumpt_DDA")

#++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++
# Initial C flux in DDA
#++++++++++++++++++++++++++++++++++++++++
tot = Fpho
a = Mu*QcD*(1+E)/tot*100
b = (Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC-FphoV)/tot*100
c = FphoV/tot*100

values = [nanmean(a),nanmean(b),nanmean(c)]
#-----------------------------------------------------------------------
# Checking (if these are not equal, that means there is Mu dependence)
#-----------------------------------------------------------------------
    
if round(a[-1],10) != round(nanmean(a),10):
    print("no good")
    Print(a)
if round(b[-1],10) != round(nanmean(b),10):
    print("no good")
    Print(b)
if round(c[-1],10) != round(nanmean(c),10):
    print("no good")
    Print(c)
#-----------------------------

figure(25)
pie(values,labels=['Dia','Dia-Chain','Chain'],autopct='%1.0f%%',startangle=90,counterclock=False,textprops={'size':23},wedgeprops=dict(edgecolor='k',linewidth=1.5))#,textprops={'size':30})
title("Inicial C flux in DDA")
sf("Pie_C_flux_initial_DDA")


print(Fpho/Fn2fixN/14*12)
print(Qn/QnC)
#++++++++++++++++++++++++++++++++++++++

if PanelPlot == 1:
    Savefig3('02\\04 DDA','LargePlot',300)

show()
