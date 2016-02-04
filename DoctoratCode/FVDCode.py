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
if __name__=='__main__':
    import matplotlib.pylab as plt
    plt.ion()
else:
    from matplotlib.pylab import ion
    ion()

fig1 = plt.figure(1)
ax1 = plt.subplot(3,1,1)
ax2 = plt.subplot(3,1,2)
ax3 = plt.subplot(3,1,3)

def eps(ntu,mu):
    return 1.-np.exp(1./mu*ntu**0.22*(np.exp(-mu*ntu**0.78)-1))

def kA(alphaa,Aa,effa,alphaw,Aw,effw=1.):
    return 1./(1./(alphaa*Aa*effa)+1./(alphaw*Aw*effw))

Qm = 68.3                # kW

ma = 2.107               # kg/s
va = 8.71                # m/s
Tina = 18.15             # C
Touta = 50.42            # C
alphaa = 94.6            # W/m^2/K
effa=1.0
rhoa = 1.0/HAPropsSI('V','T',(Tina+Touta)/2.+273.15,'P',1.013e5,'R',0.4)        #kg/m3
cpa = HAPropsSI('C','T',(Tina+Touta)/2.+273.15,'P',1.013e5,'R',0.4)             #J/kg/k
Qa = ma*cpa*(Touta-Tina) #W
print("Qa=%.3f [W]"%Qa)

vdotw =264.0                #l/min 
vw = 0.8202              #m/s
Tinw = 90.3              #C
Toutw = 86.5             #C
alphaw = 8221.6          #W/m^2/K
rhow = PropsSI('D','T',(Tinw+Toutw)/2.0+273.15,'P',2.0e5,'Water')               #kg/m3
cpw =PropsSI('C','T',(Tinw+Toutw)/2.+273.15,'P',2.0e5,'Water')                  #J/kg/K
mw = vdotw/60*rhow/1000. #kg/s               
Qw = mw*cpw*(Tinw-Toutw)
print("Qw=%.3f [W]"%Qw) 

#raw_input("press enter to continue")


cp_al=910.               #J/kg/K
rho_al=2712.0            #kg/m**3

L =0.830                 #m
G =0.085                 #m
ha2 = 0.008/2.0          #m  -jumatate din inaltimea canalului de aer
hw2 = 0.0022/2.0         #m  -jumatate din inalitmea canalului de apa

#volumul total al aerului din schimbator este egal cu:
Vola = L*G*(2*ha2)*(35+1)     #m^3  (35 sunt numarul de canale de apa)
Ata = 14.59                 #m^2
Volw = L*G*(2*hw2)*(35)       #m^3
Atw = 5.3998                #m^2
Vol_al = L*G*0.5e-3*(35+1)  #m^3 este un pic mai mare volumul aluminiului  


Ntu = kA(alphaa,Ata,effa,alphaw,Atw)/(ma*cpa) # W/W
mu = ma*cpa/(vdotw/60.*rhow/1000.*cpw) 
assert(mu<1.)

Qm = ma*cpa*(Tinw-Tina)*eps(Ntu,mu)
print("Qm=%.3f [W]"%Qm)
T1out_a = Tina + Qm/(ma*cpa)
print("Touta=%.2f [C]"%T1out_a)
T1out_w = Tinw - Qm/(vdotw/60.*rhow/1000.*cpw)
print("Toutw=%.2f [C]"%T1out_w)

#raw_input("press enter to continue")

alpha1a=alphaa*Ata
alpha1w=alphaw*Atw

nx = 100
dx = L/nx
ny = 20
dy = G/ny

mesh2d = Grid2D(dx=dx,dy=dy,nx=nx,ny=ny)
Tc = CellVariable(mesh=mesh2d,name="$T_{aer}$",value=Tina,hasOld=1)
Th = CellVariable(mesh=mesh2d,name='$T_{apa}$',value=Tina,hasOld=1)
Ts =CellVariable(mesh=mesh2d,name='$T_{perete}$',value=Tina,hasOld=1) #

Tc.constrain(Tina,mesh2d.facesBottom)
Tc.faceGrad.constrain([0],mesh2d.facesTop)

Th.constrain(Tinw,mesh2d.facesLeft)
Th.faceGrad.constrain([0],mesh2d.facesRight)

eqcold=rhoa*cpa*Vola*TransientTerm(1.,var=Tc) +ma*cpa*G*PowerLawConvectionTerm((0.0,1.),var=Tc) ==alpha1a*effa*(Ts-ImplicitSourceTerm(1.,var=Tc))
eqwal=TransientTerm(cp_al*rho_al*Vol_al,var=Ts) == alpha1w*(Th-ImplicitSourceTerm(1.,var=Ts)) + alpha1a*effa*(Tc-ImplicitSourceTerm(1.,var=Ts))
eqhot=rhow*cpw*Volw*TransientTerm(1.,var=Th) + mw*cpw*L*PowerLawConvectionTerm((1.,0.0),var=Th) == alpha1w*(Ts-ImplicitSourceTerm(1.,var=Th))

aspectratio=.5
if __name__=='__main__':
    vTc = Matplotlib2DViewer(Tc,'a. Temperatura aerului',axes=ax1,figaspect=aspectratio,datamin=Tina-0.5,datamax=Touta+0.5)
    #ax1.set_ylabel('Grosimea [m]')
    ax1.yaxis.set_ticks([0.0,0.02,0.04,0.06,0.08])
    vTs = Matplotlib2DViewer(Ts,'b. Temperatura peretelui',axes=ax2,figaspect=aspectratio,datamin=Toutw-2.,datamax=Tinw+0.5)
    ax2.set_ylabel('Grosimea [m]')
    ax2.yaxis.set_ticks([0.0,0.02,0.04,0.06,0.08])
    vTh = Matplotlib2DViewer(Th,'c. Temperatura apei',axes=ax3,figaspect=aspectratio,datamin=Toutw-2.,datamax=Tinw+0.5)
    #ax3.set_ylabel('Grosimea [m]')
    ax3.set_xlabel('Lungimea [m]')
    ax3.yaxis.set_ticks([0.0,0.02,0.04,0.06,0.08])
    viewer = MultiViewer([vTc,vTs,vTh]) #,datamin=Tina-1.,datamax=Tinw+1.)

eqs=eqcold & eqwal & eqhot
dt=.05
steps=70

Tmouta=np.zeros(steps)
Tmoutw=np.zeros(steps)
Tmsurf=np.zeros(steps)

print('begin calc...')
for step in range(steps):
    Tc.updateOld()
    Ts.updateOld()
    Th.updateOld()
    eqs.solve(dt=dt)
    Tmouta[step]=np.mean(Tc.faceValue[mesh2d.facesTop.value].value)
    Tmoutw[step]=np.mean(Th.faceValue[mesh2d.facesRight.value].value)
    Tmsurf[step]=np.mean(Ts.value)
    ax1.set_aspect(2)
    ax2.set_aspect(2)
    ax3.set_aspect(2)
    viewer.plot()
    if step%8==0:
        fig1.suptitle('timp = %.2f'%(step*dt))
        fig1.savefig('heatmap%.1f.png'%(step*dt),dpi=300.0)
    print('timp = %.2f [s]'%(step*dt))

print('final t=%.2f [s]'%(steps*dt))

T2outa=np.mean(Tc.faceValue[mesh2d.facesTop.value].value)
T2outw=np.mean(Th.faceValue[mesh2d.facesRight.value].value)
print("Tmouta=%.3f [C]"%T2outa)
print("Tmoutw=%.3f [C]"%T2outw)
fig1.suptitle('timp = %.2f [s]'%(steps*dt))
fig1.savefig('heatmap.png',dpi=300.0)

fig2 = plt.figure(2)
ax = plt.subplot(311)
time=np.array(range(steps))*dt
plt.plot(time,Tmouta,'b-')
plt.title('a. Temperatura medie de iesire a aerului')
plt.ylabel('$To_c [\degree C]$')
ax.grid(True)


ax = plt.subplot(312,sharex=ax)
plt.plot(time,Tmsurf,'g-')
plt.title('b. Temperatura medie a peretelui')
plt.ylabel('$Ts [\degree C]$')
ax.grid(True)

ax = plt.subplot(313,sharex=ax)
plt.plot(time,Tmoutw,'r-')
plt.title('c. Temperatura meide de iesire a apei')
plt.ylabel('$To_h [\degree C]$')
plt.xlabel('Timp [s]')
ax.grid(True)

plt.tight_layout()

fig2.savefig('Tmout.png',dpi=300.0)

print("gata")
