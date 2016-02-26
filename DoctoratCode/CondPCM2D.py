"""
Calculam distributia de temperatura datorita conductiei intre perete si pcm:

calculul se va face pe 1/4 din domeniul simulat considerand simetriile problemei
vezi descrierea problemei in CondPCM2D.svg

proprietatile alum 3003:
[1] http://www.makeitfrom.com/material-properties/3003-H14-Aluminum/

proprietatile pcm A25H:
[2] ..\Articole ISI\ideea placilor de racire cu PCM
"""


from fipy import *
import numpy as np
if __name__=='__main__':
    import matplotlib.pyplot as plt
    plt.ion



rhoal  = 2730.            #kg/m^3
cpal   = 890.             #J/kg.K
lal    = 160.             #W/m.K
Dal    =lal/(rhoal*cpal)  #coeficientul de difuzie

rhopcm = 810.0            #kg/m^3
lpcm   = 0.180            #W/m.K
llpcm  = 226e3            #W/kg
cppcm  =2150.0            #J/kg.K
Dpcm   =lpcm/(rhopcm*cppcm) #coef de difuzie

alpha =100.               #W/m^2.K


l=2.5e-3
h=4.2e-3
gf=0.2e-3
gpd=0.5e-3
A=l/2.*1.6e-3

m = Grid2D(nx=30,ny=30,Lx=l/2.,Ly=h/2.)
#m1d = Grid1D(nx=30,Lx=h/2.)
Tinf=15.

T=CellVariable(mesh=m,value=5.,hasOld=True)
#T1d=CellVariable(mesh=m1d,value=5.,hasOld=True)
#D=FaceVariable(mesh=m,value=Dpcm)
#X,Y = m.faceCenters
#al_mask = ((X<=gf/2.) | (Y>=h/2.-gpd))
#D.setValue(Dal,where=al_mask)
#T.constrain(15.,m.facesTop)

D=(Dal-Dpcm)*((m.x<=gf/2.)|(m.y>=(h/2.-gpd)))+Dpcm #*((m.x>gf/2.)&(m.y<(h/2.-gpd)))
DVar=CellVariable(mesh=m,name='D',value=D)

T.faceGrad.constrain(alpha/lal*(Tinf-T.faceValue)*m.faceNormals,m.facesTop)
T.faceGrad.constrain(((0.,),(0.,)),where=m.facesLeft | m.facesBottom | m.facesRight)

#T1d.faceGrad.constrain(alpha/lpcm*(Tinf-T1d.faceValue)*m1d.faceNormals,m1d.facesRight)
#T1d.faceGrad.constrain(0.*m1d.faceNormals,m1d.facesLeft)


if __name__=='__main__':
    Tview=Viewer(vars=T,title='dif coef',datamin=5.,datamax=16.)
    Dview=Viewer(vars=DVar,title='diff')
    Dview.plot()
    Dview.plotMesh('mesh2D')
    #T1dview=Viewer(vars=T1d,title='1dim ')

eq=TransientTerm()==DiffusionTerm(coeff=D)

#eq1d=TransientTerm()==DiffusionTerm(coeff=Dpcm)

dt=0.8
for i in range(200):
    T.updateOld()
    T.faceGrad.constrain(alpha/lal*(Tinf-T.faceValue)*m.faceNormals,m.facesTop)
    #T1d.updateOld()
    #T1d.faceGrad.constrain(alpha/lpcm*(Tinf-T1d.faceValue)*m1d.faceNormals,m1d.facesRight)
    eq.solve(T,dt=dt)
    #eq1d.solve(T1d,dt=dt)
    if __name__=='__main__':
        Tview.setLimits(datamin=np.min(T.value)-0.5,datamax=np.max(T.value)+0.5)
        Tview.plot()
        #T1dview.plot()



if __name__=='__main__':
    #Dview=Viewer(vars=T,title='dif coef',)
    #Dview.plot()
    raw_input('enter..')

pass

               
