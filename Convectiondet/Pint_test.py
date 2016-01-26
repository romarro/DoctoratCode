# -*- coding: utf-8 -*-
"""
Created on Tue Jun 03 19:41:06 2014

@author: Vlad
"""
from __init__ import um,Q_
import CoolProp.CoolProp as fld
import uncertainties as un
import numpy as np


d=Q_(158,um.liters/um.minutes)
cp=Q_(4.21,um.kiloJ/(um.kilogram*um.K))
dens=Q_(1000,um.kg/(um.m**3))
Tin=Q_(60,um.degC).to(um.degK)
Tout=Q_(58,um.degC).to(um.degK)
Q=d*dens*cp*(Tin-Tout)

@un.wrap
def dens(tmp,press):
    return fld.Props('D','T',tmp,'P',press,'Water')
    
    
inw_dens=np.vectorize(dens)
inw_um_dens=um.wraps(um.kg/um.m**3,[um.degK,um.kPa],strict=False)(inw_dens)

class tst_fld():
    
    def __init__(self,T=0*um.degC,p=1*um.atm,Rh=0):
        if(not isinstance(p,Q_)):
            self.__p=Q_(p,um.atm)
        else:
            self.__p=p
        self.__T=T               
   
    @um.wraps(um['kg/m**3'],[None,um.degK],False)
    @np.vectorize    
    @un.wrap
    def dens(self,tmp):
        p=(self.__p.to(um.kPa)).magnitude
        return fld.Props('D','T',tmp,'P',p,'Water')
   

print 'deb=',d
print 'cp=',cp
print 'dens',dens
print 'Tin=',Tin.to('degC'),' Tout=',Tout.to('degC')
print 'Q=',Q.to(um.kW)
print 'Tin-Tout=',(Tin-Tout)
print inw_um_dens([un.ufloat(20,0.1),un.ufloat(30,0.2),un.ufloat(40,0.01)]*um.degC,1*um.atm)
print inw_um_dens(np.array([20+273.15,30+273.15,40+273.15]),100)

print '--------------test ums with un in __init__-------------'
water=tst_fld()
print (water.dens([20,30,40]*um.degC)*(158*um['l/min'])).ito_base_units()
print water.dens([un.ufloat(20,0.1),un.ufloat(30,0.2),un.ufloat(40,0.01)]*um.degC)
print water.dens([20+273.15,30+273.15,40+273.15])
print '--------------test ums with un in __init__-------------'
water1=tst_fld(20.,1,0.5)