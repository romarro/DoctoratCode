# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 08:06:34 2014

@author: Vlad
"""
import uncertainties as un

def uequal(x,y):
    def ininterval(a,minimum,maximum):
        if a<=maximum and a>=minimum:
            return True
        else:
            return False
            
    if isinstance(x,un.UFloat) and isinstance(y,un.UFloat):
        if x.n>y.n:
            swap=x
            x=y
            y=swap
        
        return ininterval(x.n+x.s,y.n-y.s,y.n+y.s)
        
    if isinstance(x,float) and isinstance(y,un.UFloat):
        return ininterval(x,y.n-y.s,y.n+y.s)
    if isinstance(x,un.UFloat) and isinstance(y,float):
        return ininterval(y,x.n-x.s,x.n+x.s)



if __name__=='__main__':
    y=un.ufloat(1.0,0.5)
    x=un.ufloat(1.4,0.2)
    print 'TRUE ASSERT : {:P} == {:P} este {:}'.format(x,y,uequal(x,y))
    x=un.ufloat(3,0.5)
    y=un.ufloat(2,0.2)
    print 'FALSE ASSERT: {:P} == {:P} este {:}'.format(x,y,uequal(x,y))
        
        
            
        