'''
Created on Apr 16, 2019
From D002 01 00
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

def comp(Nh):
    #-------------------------------
    # Light energy coverage (LEC)
    #-------------------------------
    LEC = arange(0,1+0.01,0.01)
    #-------------------------------
    
    CN = 6.6
    Nv = 4 #Number of the vegetative cells per heterocysts
    #Nh = 2 #Number of the heterocysts per host diatom
    
    R = arange(6,11+0.1,0.1)
    
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
    
    QnV = QcV/CN #(pmol N cell-1)
    QnH = QcH/CN #(pmol N cell-1)
    QnD = QcD/CN #(pmol N cell-1)
    
    QnC = QnV + QnH
    Qn = QnV + QnH + QnD
    
    #Maybe I can make it growth rate array and use it for the x-axis.
    Mu = 0.51 #(d-1) Growth rate (Follett 2018)
    
    Fn2fixN = Mu*(QnV + QnH + QnD) #(pmol N cell-1 d-1)
    
    E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate
    
    YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
    #-----------------------------
    # Using LEC
    #-----------------------------
    YresN2fix = YresN2fix*(1-LEC)
    #-----------------------------
    
    Fn2fixC = Fn2fixN*YcnN2fix + Fn2fixN*YresN2fix #(pmol C cell-1 d-1) C consumption for N2 fixation
    
    Fpho = Mu*(QcV + QcH + QcD)*(1+E) + Fn2fixC #(pmol C cell-1 d-1) Photosynthesis rate
    
    FphoV = Fpho*QnV/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by vegetative cells
    FphoD = Fpho*QnD/(QnV + QnD) #(pmol C cell-1 d-1) Photosynthesis rate by diatom cells
    
    PanelPlot = 0
    if PanelPlot == 1:
        figure(1,figsize=(25,12))
        rcParams.update({'font.size': 15})
        
    def sup(n):
        if PanelPlot == 1:
            subplot(2,4,n)
        else:
            figure(n)
        xlabel('LEC$_{Het}$ ($\%$)')
    
    def sf(name):
        if PanelPlot != 1:
            Savefig3('02\\04 DDA',name,300)
    
    ColorN = ('','#FFCC00','#DF9900','#C06600','#A03300','#800000') #From the Sheet 1 in "C:\Users\keiin\Google Drive\0 Paper\71 Modeling DDA\Excel\07 ColorForPlot.xlsx"
    ColorC = ('','#FFCC00','#DF9900','#C06600','#A03300','#800000')
    ColorSC = ('','#FFCC00','#DF9900','#C06600','#A03300','#800000')
    
    sup(12)
    FcCostChain = Mu*QcV+Mu*QcV*E+Mu*QcH+Mu*QcH*E+Fn2fixC  #C cost by chain including N supply to diatom
    Output = (FcCostChain-FphoV)/FcCostChain*100
    if Nh == 1:
        plot([50,50],[0,100],color="white") #To make the y range same Fig 5B
    plot(LEC*100,Output,color=ColorC[Nh],label='Chain x '+str(Nh))
    title('C supply from Diatom',y=1.02)
    ylabel('C supply ($\%$)')
    if Nh == 5:
        sf('LED_Infuence_on_C_supply_from_Diatom')
    
    print(Output[0]-Output[-1])
    
comp(1)
comp(2)
comp(3)
comp(4)
comp(5)


show()
