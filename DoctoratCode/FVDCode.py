"""
Pentru a vedea cum variaza temperatura in lungul aripioarei 
ondulate vom face un calcul 2D pe un racitor si a vedea daca 
presupunerea noastra ca temperatura este constanta se poate 
aplica. In plus trebuie sa vedem si care sunt diferentele de 
temperatura pentru aer si apa

Ne vom folosi de date experimentale date de RA 28788-0
in special pozitia 23 pentru testarea  metodei de calcul
"""
from __future__ import division
import numpy as np
from CoolProp.CoolProp import HAPropsSI,PropsSI
from fipy import *


Qm = 68.3                # kW
ma = 2.107               # kg/s
va = 8.71                # m/s
Tina = 18.15             # C
Touta = 50.42            # C
alphaa = 94.6            # W/m^2/K

mw =264.0                #l/min
vw = 0.8202              #m/s
Tinw = 90.3              #C
Toutw = 86.5             #C
alphaw = 8221.6          #W/m^2/K

cp_al=910.               #J/kg/K
rho_al=2712.0            #kg/m**3

L =0.830                 #m
G =0.085                 #m
ha2 = 0.008/2.0          #m  -jumatate din inaltimea canalului de aer
hw2 = 0.0022/2.0         #m  -jumatate din inalitmea canalului de apa

#volumul total al aerului din schimbator este egal cu:
Vola=L*G*(ha2)*(35+1)   #m^3  (35 sunt numarul de canale de apa)
Ata=14.59                 #m^22
Volw=L*G*(hw2)*(35)     #M^3
Atw=5.3998                #m^2

rhoa = 1.0/HAPropsSI('V','T',(Tina+Touta)/2.+273.15,'P',1.013e5,'R',0.4)        #kg/m3
cpa = HAPropsSI('C','T',(Tina+Touta)/2.+273.15,'P',1.013e5,'R',0.4)             #J/kg/k

rhow = PropsSI('D','T',(Tinw+Toutw)/2.0+273.15,'P',2.0e5,'Water')               #kg/m3
cpw =PropsSI('C','T',(Tinw+Toutw)/2.+273.15,'P',2.0e5,'Water')                  #J/kg/K

print "debit aer recalc"

alpha1a=alphaa*Ata/Vola
alpha1w=alphaw*Atw/Volw
effa=0.98

nx = 200
dx = L/nx
ny = 20
dy = G/ny

mesh2d = Grid2D(dx=dx,dy=dy,nx=nx,ny=ny)
Tc = CellVariable(mesh=mesh2d,name="$T_{aer}$",value=Tina,hasOld=1)
Th = CellVariable(mesh=mesh2d,name='$T_{apa}$',value=Tinw,hasOld=1)
#Ts = CellVariable(mesh=mesh2d,name="$T_{perete}$",value=80.0) #C constant in prima faza
Ts =CellVariable(mesh=mesh2d,name='$T_{perete}$',value=Tina,hasOld=1) #

Tc.constrain(Tina,mesh2d.facesBottom)
Tc.faceGrad.constrain([0],mesh2d.facesLeft)
Tc.faceGrad.constrain([0],mesh2d.facesRight)
Tc.faceGrad.constrain([0],mesh2d.facesTop)

Th.constrain(Tinw,mesh2d.facesLeft)
Th.faceGrad.constrain([0],mesh2d.facesBottom)
Th.faceGrad.constrain([0],mesh2d.facesTop)
Th.faceGrad.constrain([0],mesh2d.facesRight)

Ts.faceGrad.constrain([0],mesh2d.facesLeft)
Ts.faceGrad.constrain([0],mesh2d.facesBottom)
Ts.faceGrad.constrain([0],mesh2d.facesTop)
Ts.faceGrad.constrain([0],mesh2d.facesRight)

eqcold=rhoa*cpa*TransientTerm(1.,var=Tc) + rhoa*cpa*PowerLawConvectionTerm((0.0,va),var=Tc) ==alpha1a*effa*(Ts-ImplicitSourceTerm(1.,var=Tc))
eqwal=TransientTerm(cp_al*rho_al,var=Ts) == alpha1w*(Th-ImplicitSourceTerm(1.,var=Ts)) + alpha1a*effa*(Tc-ImplicitSourceTerm(1.,var=Ts))
eqhot=rhow*cpw*TransientTerm(1.,var=Th) + rhow*cpw*PowerLawConvectionTerm((vw,0.0),var=Th) == alpha1w*(Ts-ImplicitSourceTerm(1.,var=Th))

if __name__=='__main__':
    #viewer = MultiViewer([Viewer(Tc),Viewer(Th),Viewer(Ts)])
    viewer = Viewer([Tc,Th,Ts],datamin=Tina-1.,datamax=Tinw+1.)

eqs=eqcold & eqwal & eqhot
dt=.4
for step in range(20):
    eqs.solve(dt=dt)
    viewer.plot()
    raw_input("time %s next step:"%(step*dt))

#eqs.solve(dt=0.5)
#eqcold.solve(dt=1.)
#for step in range (10):
#    eqcold.solve(dt=.3)
#    eqhot.solve(dt=.3)
#    if __name__=='__main__':
#        viewer.plot()
#        raw_input("next step")

raw_input("press a key")

