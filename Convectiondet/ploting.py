# -*- coding: utf-8 -*-
"""
Created on Fri Jun 06 08:03:39 2014

@author: Vlad


"""
import pint as pint
from __init__ import um,Q_
import xlrd as xls
import matplotlib.pyplot as plt
import uncertainties as un
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as st
from itertools import product

def unitsPrnt(x):
    y=x.units.__class__(dict((x._REGISTRY.get_symbol(key), value) for key, value in x.units.items()))
    return y

def uplot(x,y,fmt='.',func=None,sg=2,**kwargs):
    """
    prototipul functiei::
        def func(x,*args):
            return ...
    
    """
    noerr=kwargs.pop('noerr',False)
    
    
    def anom(arr):
        return np.array([k.nominal_value for k in arr])
    def astd(arr):
        return np.array([k.std_dev for k in arr])    

    
    if isinstance(x,Q_):
        x=x.tolist()
        xn=np.array([s.magnitude for s in x])
    else:
        xn=x
    if isinstance(y,Q_):
        y=y.tolist()
        yn=np.array([s.magnitude for s in y])
    else:
        yn=y
    
    #in errorbar se introduc erorile relative si nu cele absolute
    xv=anom(xn)
    xs=astd(xn)
    yv=anom(yn)
    ys=astd(yn)
    
    slim=kwargs.pop('slim',max(xv))
    ilim=kwargs.pop('ilim',min(xv))
    
    pl=plt.errorbar(xv,yv,ys,xs,fmt,ecolor='black',**kwargs)
          
    if func!=None:
        #determinam coeficientii care fiteaza curba cu functia data
        popt,pcov=curve_fit(func,xv,yv,sigma=ys/yv)
        
        #determinam erorile coeficientilor functiei
        sigmas=np.sqrt([pcov[i,i] for i in range(len(pcov))])
        maps=map(np.array,product([1,2],repeat=len(sigmas)))
        p=map(lambda x:np.power(-1,x),maps)
        scoefs=map(tuple,popt+p*sigmas)
        xk=np.linspace(ilim,slim,100)        
        vals=[func(xk,*scoef) for scoef in scoefs]
        fstd=np.std(vals,axis=0)              
        plt.plot(xk,func(xk,*tuple(popt)),'g-',label='aprox.')
        if not noerr:
            plt.plot(xk,func(xk,*tuple(popt))+sg*fstd,'r:',
                     label=r'+{:}$\cdot\sigma$'.format(sg))
            plt.plot(xk,func(xk,*tuple(popt))-sg*fstd,'r:',
                     label=r'-{:}$\cdot\sigma$'.format(sg))
        return (pl,[un.ufloat(x[0],x[1]) for x in zip(popt,np.sqrt(np.diag(pcov)))])
    return pl
    



#-------------------pentru testare-------------------------------
def gendatalog(lim=(10,100),pcts=10,smaxx=0.5,smaxy=0.3):
    np.random.seed(0)    
    x=np.linspace(lim[0],lim[1],pcts)
    y=2.45*np.log(x)
    sigmay=np.abs(np.random.normal(0,np.ones(len(y))*smaxy))
    np.random.seed(0)    
    y=np.random.normal(y,sigmay)
    sigmax=np.abs(np.random.normal(0,np.ones(len(x))*smaxx))
    x=np.random.normal(x,sigmax)
    return (x,y,sigmax,sigmay)
 


def loadDta(filename,offset=(1,0),rows=20,cols=13,index=0):
    wb=xls.open_workbook(filename)
    sh=wb.sheet_by_index(index)
    data=np.array([[sh.cell_value(r,c) for c in range(offset[1],cols+offset[1])] for r in range(offset[0],rows+offset[0])])
    data[data=='']='0'
    data.astype(float)    
    return data   

if __name__=='__main__':
    x=13.0*um.parse_expression('W/m/K')
    x,y,sigmax,sigmay=gendatalog(pcts=20)
    ux=[un.ufloat(v,er) for v,er in zip(x,sigmax)]*um.l/um.s
    uy=[un.ufloat(v,er) for v,er in zip(y,sigmay)]*um.kW
    #plt.plot(x,y,'b.')
    #plt.errorbar(x,y,sigmax,sigmay,fmt='o')
    pl,coefs=uplot(ux,uy,func=lambda x,a,b:a*x**b,fmt='k.',noerr=True)
    plt.ylabel('Flux termic [{:~}]'.format(unitsPrnt(uy)))
    plt.xlabel('debit aer [{:~}]'.format(unitsPrnt(ux)))
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.show()
    print coefs