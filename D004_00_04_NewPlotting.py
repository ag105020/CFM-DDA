'''
Created on Apr 16, 2019
From D000 00 09 
Computing how many vegetative cells needed
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

#print(QcH/QcC)

QnV = QcV/CN #(pmol N cell-1)
QnH = QcH/CN #(pmol N cell-1)
QnD = QcD/CN #(pmol N cell-1)

QnC = QnV + QnH
Qn = QnV + QnH + QnD

#print(QcV/QnV,QcH/QnH,QcD/QnD)

#Maybe I can make it growth rate array and use it for the x-axis.
MuMax = 0.31 #(d-1) Maximum growth rate of Richelia; Villareal 1990 (from Foster 2011)
MuStep = 0.01
Mu = arange(0.0,MuMax+MuStep,MuStep) #(d-1) Growth rate
Mu = 0.31 #(d-1) Maximum growth rate of Richelia; Villareal 1990 (from Foster 2011)

Fn2fixN = Mu*(QnV + QnH + QnD) #(pmol N cell-1 d-1) 

E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate

YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron

#print(YcnN2fix,YresN2fix)

Fn2fixC = Fn2fixN*YcnN2fix + Fn2fixN*YresN2fix #(pmol C cell-1 d-1) C consumption for N2 fixation

Fpho = Mu*(QcV + QcH + QcD)*(1+E) + Fn2fixC #(pmol C cell-1 d-1) Photosynthesis rate

FphoV = Fpho*QnV/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by vegetative cells
FphoD = Fpho*QnD/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by diatom cells


##########################
# Computation of ratios
##########################
#Mu = arange(0.0,MuMax+MuStep,MuStep) #(d-1) Growth rate
Rmu = 1

CsVeg = FphoV/(Nv*Nh)   #(pmol C cell-2 d-1) C supply per vegetative cells
CdVeg = Mu*QcV*(1+E)/(Nv*Nh)   #(pmol C cell-1 d-1) Growth related C demand per vegetative cells 
CdNVeg = Fn2fixC/Qn*QnV/(Nv*Nh)*Rmu   #(pmol C cell-1 d-1) N fixation C cost per vegetative cells
CdHet = Mu*QcH*(1+E)/Nh   #(pmol C cell-1 d-1) Growth relatd C demand per heterocyst
CdNHet = Fn2fixC/Qn*QnH/Nh*Rmu  #(pmol C cell-1 d-1) N fixation C cost per heterocyst

CdNDia = Fn2fixC/Qn*QnD*Rmu

NVeg0 = (CdHet + CdNHet)/(CsVeg - CdVeg - CdNVeg)
print(NVeg0)
NVeg2 = (CdHet + CdNHet + CdNDia/Nh)/(CsVeg - CdVeg - CdNVeg)
print(NVeg2)

nVeg = arange(0,36+0.1,0.1)
Demand1 = (CdVeg+CdNVeg)*nVeg*Nh + (CdHet + CdNHet)*Nh
Demand2 = (CdVeg+CdNVeg)*nVeg*Nh + CdNDia + (CdHet + CdNHet)*Nh
Supply = CsVeg*nVeg*Nh

####################
# Plotting
####################

def sf(name):
        Savefig3('02\\04 DDA',name,300)


###################################
nVeg = arange(0,104+0.1,1)
Demand1 = (CdVeg+CdNVeg)*nVeg*Nh + (CdHet + CdNHet)*Nh
Demand2 = (CdVeg+CdNVeg)*nVeg*Nh + CdNDia + (CdHet + CdNHet)*Nh
Supply = CsVeg*nVeg*Nh


figure(5)
plot(nVeg,(nVeg*Nh*Vv0+Vh)/Vd*100)
plot(nVeg,ones(size(nVeg))*100,':',color='k')
title('Space occupation by the trichomes', y=1.02)
xlabel('Vegetative cells / Heterocyst')
ylabel('Space occupation (%)')
sf('Space_occupation')


figure(12)
plot(nVeg,Supply - Demand2,color="#2CA02C")
plot(nVeg,zeros(size(nVeg)),':',color='k')
#plot(nVeg,(Supply/Demand1)*100,color="#2CA02C")
#plot(nVeg,ones(size(nVeg))*100,':',color='k')
title('Difference in C supply and demand', y=1.02)
xlabel('Vegetative cells / Heterocyst')
ylabel('Supply $-$ Demand (pmol C d$^{-1}$)')
sf('Supply-Demand')

a = Supply/Demand2

show()
