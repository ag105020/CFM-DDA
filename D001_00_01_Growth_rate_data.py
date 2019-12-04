'''
Created on Apr 16, 2019
Beginning DDA Modeling
Kei215-4~7
@author: keiin
'''

from pylab import *
from FigSetting2 import *
from Savefig3 import *

def sf(name):
        Savefig3('02\\04 DDA',name,300)

Data = array([0.315769231,0.311,0.5118])   #From Follett 2018 ("05 Growth of different diazotrophs.xlsx")
SD = array([0.176279207,0.127204372,0.149]) #Same as "Data"
Names = ['Tricho.','Croco.','DDA']

bar(Names,Data,yerr=SD,color='#2CA02C',edgecolor='k',linewidth=2,error_kw=dict(lw=2, capsize=12, capthick=2))#width=0.2,edgecolor='k',linewidth=1.5)
title('Measured growth rate',y=1.02)
ylabel('$\mathit{\mu}$ (d$^{-1}$)')
ylim(top=0.7)
sf("Growth rate data")

show()
