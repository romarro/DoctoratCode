# -*- coding: utf-8 -*-
"""
Created on Fri Jun 06 08:03:39 2014

@author: Vlad


"""
import xlrd as xls
import matplotlib.pyplot as plt
import uncertainties as un
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as st

def pdense(x, y, sigma, M=1000):
    """ 
    Plot probability density of y with known stddev sigma.
    Functia este luata din `visualising Uncertainty 
    <http://blog.enthought.com/general/visualizing-uncertainty/#.U5FK4PnCaSq>`_   
    """
    assert len(x) == len(y) and len(x) == len(sigma)
    N = len(x)
    # TODO: better y ranging
    ymin, ymax = min(y - 2 * sigma), max(y + 2 * sigma)
    yy = np.linspace(ymin, ymax, M)
    a = [np.exp(-((Y - yy) / s) ** 2) / s for Y, s in zip(y, sigma)]
    A = np.array(a)
    A = A.reshape(N, M)
    plt.imshow(-A.T, cmap='gray', aspect='auto',
               origin='lower', extent=(min(x)[0], max(x)[0], ymin, ymax))
    plt.title('Density plot')
 

def loadDta(filename,offset=(1,0),rows=20,cols=13,index=0):
    wb=xls.open_workbook(filename)
    sh=wb.sheet_by_index(index)
    data=np.array([[sh.cell_value(r,c) for c in range(offset[1],cols+offset[1])] for r in range(offset[0],rows+offset[0])])
    data[data=='']='0'
    data.astype(float)    
    return data   

if __name__=='__main__':
    
    #data=loadDta('UncertaintiesPlots_test.xlsx',offset=(1,0),rows=20,cols=4)
    def f(x):
        return 2.45*np.log(x)
    
    #fitting curve
    def ft(x,a):
        return a*np.log(x)
        
    x=np.linspace(10,30,1000)
    y2=f(x)
    np.random.seed(0)
    y=np.random.normal(y2,0.1)
    sigma=np.abs(np.random.normal(np.zeros(len(y)),0.025))
    coefs,cov=curve_fit(ft,x,y,sigma=sigma)
    err_y=y-ft(x,coefs[0])
    s_err=np.sum(np.power(err_y,2))
    mean_x=np.mean(x)
    sumxx=((x-mean_x)**2).sum()
    n=x.size
    confid=st.t.ppf(1-0.025,n-2)*np.sqrt(s_err/(n-2))*np.sqrt(1./n+(x-mean_x)**2/sumxx)    
    
    x2=np.atleast_2d(x).T

    pdense(x2,y,confid)
    plt.plot(x,y,'k^')
    plt.plot(x,y2,'b-')
    plt.plot(x,y-abs(confid),'r--')
    plt.plot(x,y+abs(confid),'r--')