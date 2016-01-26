# -*- coding: utf-8 -*-
"""
Created on Fri Aug 08 10:14:14 2014

@author: Vlad
"""
import numpy as np
from scipy.optimize import curve_fit
import pylab as plt

def func((x,y),a,b,c):
    return (a*(x**b)*(y**c)).ravel()
    
x=np.linspace(3000,6000,10)
y=np.linspace(0.78,1.19,5)
x,y=np.meshgrid(x,y)

f=func((x,y),0.125,-0.096,0.605)

fnoisy=f+0.0002*np.random.normal(size=f.shape)

popt,pcov=curve_fit(func,(x,y),fnoisy)

fig,ax=plt.subplots(1,1)
ax.hold(True)
ax.contour(x,y,fnoisy.reshape(5,10),8,color='k')
ax.contour(x,y,f.reshape(5,10),8,color='r')

"""
#-------------------------------------------------------------
#            cream sirurile pentru linear model din sklearn
#-------------------------------------------------------------
Nus.append(ras[-1].Nu.magnitude) #valoarea functiei
fh=ras[-1].h/ras[-1].p
fh.to(um.dimensionless)
Fh=fh.magnitude*np.ones(ras[-1].Nu.size)

#sirul variabilelor un array 2D [Re,Fh]
if X==None:
    X=np.array([uns.nominal_values(ras[-1].Re.magnitude),Fh]).T
else:
    tmp=np.array([uns.nominal_values(ras[-1].Re.magnitude),Fh]).T
    X=np.append(X,tmp,axis=0)
"""