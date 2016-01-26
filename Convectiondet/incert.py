# -*- coding: utf-8 -*-
"""
Created on Mon May 19 21:20:06 2014

@author: Vlad
"""
import numpy as np
import scipy.stats as st
import uncertainties as un
import pint as um
import warnings
DEBUG=False
def ucreate(data,ierr=0.,conf=95.,unit=None):
    """
    Creaza un obiect ufloat cu incertitudinea data de intervalul de confidenta
    dorit si adauga si incertitudinea instrumentului 
    
    Parametrii
    ----------
    data: list, tuble or a numpy array
        sirul de date pentru care se calculeaza incertitudinea
    ierr: float 
        eroarea instrumentului la aceeasi confidenta ca si cea actuala
        default este zero
    conf: float
        confidenta pentru care se calculeaza incertitudinea
        default=95
    unit: pint.Quantity 
        unitatea de masura pentru variabila care se creaza.
        default=None
        
    Intoarce
    -------
    out: ufloat or pint.Measurement or pint.Quantity
        intoarce un float avand si incertitudinea, sau o masuratoare folosind
        pachetul pint pentru cantitati fizice

    Observatii
    -----
    Daca data este de tip float sau int intoarce un ufloat cu eroarea data numai
    de instrument, sau daca sirul are o deviatie standard egala cu 0.0 atunci intoarce
    incertitudinea cu eroarea data de instrument.    
    In cazul in care data nu este o valoare numerica arunca o exceptie
    Daca unit este diferit de None, acesta trebuie sa fie o unitate de masura
    definita de pachetul pint si in acest caz se va intoarce un obiect pint.Measurement
    
    References
    ----------
    [1] R. Tatara and G. Lupia, “Assessing heat exchanger performance data 
    using temperature measurement uncertainty,” Int. J. Eng. Sci. …, vol. 3,
    no. 8, pp. 1–12, 2011.
    
    """    
    warnings.filterwarnings(action='error')    
    if isinstance(data,(np.ndarray,list,tuple)):
        dta=np.array(data,dtype='f')
    else:
        if isinstance(data,(float,int)):
            return un.ufloat(data,ierr)
        else:
            raise Exception('valoare numerica')
    
    xm=np.mean(dta)
    if dta.size>1:
        xstd=np.std(dta,ddof=1)
    else:
        xstd=0.0
    if xstd==0.0 or np.isnan(xstd):
        return un.ufloat(xm,ierr)
        
    xskew=st.skew(dta)
    xint=st.t.interval(conf/100.0,dta.size-1,loc=xm,scale=xstd/np.sqrt(dta.size))
    global DEBUG
    if DEBUG:
        print u'\tmean={:.3f}'.format(xm)
        print u'\tstd={:.3f}'.format(xstd)
        print u'\tskewness={:.3f}'.format(xskew)
        print u'\tstd@95%_min={:.3f}'.format(xm-xint[0])
        print u'\tstd@95%_max={:.3f}'.format(xint[1]-xm)
    xstd=xm-xint[0]
    
    try:
        return un.ufloat(xm,np.sqrt(ierr**2+xstd**2))
    except RuntimeWarning:
        print 'xm=',xm,'ierr=',ierr,'xstd=',xstd
        print 'dta:',dta

def __main():
    """
    Testing code....
    """
    global DEBUG
    DEBUG=True
    print 'running tests....'
    print 'ucreate test:'
    print '\tsample array: [1.,1,2,1.5,1.2,1.7,0.9,0.98]'    
    x=ucreate(np.array([1.,1,2,1.5,1.2,1.7,0.9,0.98]))
    print u'\tx={:P}'.format(x)
    print 'ucreate test:'
    print '\tsample list: [1.,1,2,1.5,1.2,1.7,0.9,0.98]'    
    x=ucreate([1.,1,2,1.5,1.2,1.7,0.9,0.98])
    print u'\tx={:P}'.format(x)
    print 'ucreate test:'
    print '\tsample tuple: (1.,1,2,1.5,1.2,1.7,0.9,0.98)'    
    x=ucreate((1.,1,2,1.5,1.2,1.7,0.9,0.98))
    print u'\tx={:P}'.format(x)
    print 'ucreate test:'
    print '\t data o singura valoare,ierr=0.2'
    x=ucreate(2,ierr=0.2)
    print u'\tx={:P}'.format(x)
    print 'ucreate test:'
    print '\t list [3.,3,3,3,3,3,3,3],ierr=0.2'
    x=ucreate([3.,3,3,3,3,3,3,3],ierr=0.2)
    print u'\tx={:P}'.format(x)
    
    DEBUG=False
    
if __name__ == '__main__':
    __main()

    
    
    