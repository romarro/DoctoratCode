# -*- coding: utf-8 -*-
"""
Created on Wed May 28 08:53:41 2014

@author: Vlad
"""
from __init__ import um,Q_
import CoolProp.CoolProp as fld
import numpy as np
import uncertainties as un

class Water():
    """
    implementeaza fluidul apa
    """
    def __init__(self,conc=0):
        """
        initializeaza clasa Water
        
        Parametrii
        ----------
        conc:float
            concentratia in procente
        
        Observatii
        ----------
        pentru o concentratie de 50% conc=0.5 implicit concentratia va fi de zero   
        pentru a nu produce timpi de executie mari, la initializare se creaza si 
        functii pentru calcul vectorial si cu incertitudini
        """
        self.__defpress=2.e5 #in Pa
        self.__conc=conc
        if self.__conc==0.0:
            self.__fldName='Water'
        else:
            self.__fldName=r'EG-{con}%'.format(con=str(int(self.__conc*100.0)))
            
        #intializam functiile wrap astfel incat sa se poata folosii cu ufloat
        def __w_dens(t):
            t=float(t) #obligam conversia la float chiar daca este intreg
            return fld.PropsSI('D','T',t,'P',self.__defpress,self.__fldName)
        def __w_cp(t):
            t=float(t)
            return fld.PropsSI('C','T',t,'P',self.__defpress,self.__fldName)
        def __w_cond(t):
            t=float(t)
            return fld.PropsSI('L','T',t,'P',self.__defpress,self.__fldName)
        def __w_visc(t):
            t=float(t)
            return fld.PropsSI('V','T',t,'P',self.__defpress,self.__fldName)
        
        self.__vudens=np.vectorize(un.wrap(__w_dens))   #kg/m3
        self.__vucp=np.vectorize(un.wrap(__w_cp))       #J/kg.K
        self.__vucond=np.vectorize(un.wrap(__w_cond))   #W/m.K
        self.__vuvisc=np.vectorize(un.wrap(__w_visc))   #Pa.s
    
        
    @um.wraps(um.parse_expression('kg/m**3'),[None,um.degK])
    def densit(self,temp):
        """
        calculeaza densitatea apei
        
        Parametrii
        ----------
        temp: float, uncertainties.UFloat or ndarray
            temperatura apei, in [C]
                        
        out: float,nu.UFloat ndarray
            densitatea in [kg/m3]
            
        Observatii
        ----------
        daca temperatura este ndarray va intoarce ndarray.
        daca temperatura este UFloat se va intoarce un UFloat un numar cu incertitudine
        
        Exemplu
        -------
        Ex 1:        
            import numpy as np
            water=Water(0.5)
            temps=np.array([10,20,30,40.])
            dens=water.Densit(temps)
            print 'dens=',dens
            
            la consola:
            dens=[1070.0208845 ,  1064.92866283,  1059.38820425,  1053.44076737]
        Ex 2:
            import numpy as np
            import uncertainties as un
            water=Water() #apa pura
            temps=np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)])
            dens=water.Densit(temps)
            print 'dens=',dens
            
            la consola:
            dens=[999.7496259274476+/-0.0017598602373177705 
                  998.2523477831393+/-0.004131876524661046]
        """
                   
        return self.__vudens(temp)
    
    @um.wraps(um.parse_expression('J/(kg*degK)'),[None,um.degK])    
    def cp(self,temp):
        """
        Calculeaza caldura specifica la presiune constanta
        
        Parametrii
        ----------
        temp: float, UFloat or ndarray
            temperatura in [C]
        
        out: float, UFloat or ndarray
        
        Observatii
        ----------
        daca temperatura este ndarray va intoarce ndarray.
        daca temperatura este UFloat se va intoarce un UFloat un numar cu incertitudine
        """
        return self.__vucp(temp)
   
    @um.wraps(um.parse_expression('W/m/degK'),[None,um.degK])     
    def cond(self,temp):
        
        return self.__vucond(temp)
    
    @um.wraps(um.parse_expression('Pa*s'),[None,um.degK])
    def visc(self,temp):
        
        return self.__vuvisc(temp)
    
    def Pr(self,temp):
        
        cp=self.cp(temp)
        v=self.visc(temp)
        con=self.cond(temp)
        
        return (cp*v/con).to(um.dimensionless)
        
class Air():
    """calculeaza propritetatile aerului"""
    def __init__(self,Tabs=Q_(20.,um.degC),Pabs=1*um.atm,Rh=0.0):
        """
        initializeaza proprietatile aerului in functie de umiditiatea acestuia
        
        Parametrii
        ----------
        Tabs: float, pint.Quantity or UFloat
            temperatura de referinta a mediului ambiant [C]
        Pabs: float, pint.Quantity or UFloat
            presiunea de referinta a mediului ambiant [atm]
        Rh:float, pint.Quantity or UFloat [0,1]
            umiditatea relativa a mediului ambiant
            
        Observatii
        ----------
        Se considera ca aerul va suferi transformari termodinamice fara schimbarea
        cantitatii de umiditate (fara a condensa umiditatea) astfel umiditatea
        absoluta odata calculata la initializare nu se mai modifica 
        
        """
        if not isinstance(Tabs,Q_):
            Tabs=Q_(Tabs,'degC')        
        if not isinstance(Pabs,Q_):
            Pabs=Q_(Pabs,'atm')
        if isinstance(Rh,Q_):
            Rh=Rh.magnitude
            
        self.__tabs=Tabs        
        self.__pabs=Pabs
        self.__rh=Rh
        
        @um.wraps(um.parse_expression('kg/kg'),[um.degK,um.Pa,None])
        @un.wrap
        def __w_abs(t,p,rh):
            t=float(t)
            p=float(p)
            rh=float(rh)
            return fld.HAPropsSI('W','T',t,'P',p,'R',rh)
        
        self.__abs=__w_abs(self.__tabs,
                           self.__pabs,self.__rh)
        
        def __w_dens(t,p,ab):
            t=float(t)
            p=float(p)
            ab=float(ab)
            return 1.0/fld.HAPropsSI('V','T',t,'P',p,'W',ab)
        def __w_cp(t,p,ab):
            t=float(t)
            p=float(p)
            ab=float(ab)
            return fld.HAPropsSI('C','T',t,'P',p,'W',ab)
        def __w_cond(t,p,ab):
            t=float(t)
            p=float(p)
            ab=float(ab)
            return fld.HAPropsSI('K','T',t,'P',p,'W',ab)
        def __w_visc(t,p,ab):
            t=float(t)
            p=float(p)
            ab=float(ab)
            print 't=',t
            print 'p=',p
            print 'ab=',ab
            rez=fld.HAPropsSI('M','T',t,'P',p,'W',ab)
            print rez
            return rez
        
        self.__vndens=np.vectorize(un.wrap(__w_dens))
        self.__vncp=np.vectorize(un.wrap(__w_cp))
        self.__vncond=np.vectorize(un.wrap(__w_cond))
        self.__vnvisc=np.vectorize(un.wrap(__w_visc))
    
    @um.wraps(um.parse_expression('kg/m**3'),[None,um.degK,um.Pa])            
    def densit(self,t,p=None):
        if p==None:
            p=self.__pabs.to('Pa').magnitude
        return self.__vndens(t,p,self.__abs.magnitude)
        
    @um.wraps(um.parse_expression('J/kg/degK'),[None,um.degK,um.Pa])
    def cp(self,t,p=None):
        if p==None:
            p=self.__pabs.to('Pa').magnitude
        return self.__vncp(t,p,self.__abs.magnitude)
    
    @um.wraps(um.parse_expression('W/m/degK'),[None,um.degK,um.Pa])
    def cond(self,t,p=None):
        if p==None:
            p=self.__pabs.to('Pa').magnitude
        return self.__vncond(t,p,self.__abs.magnitude)
        
    @um.wraps(um.parse_expression('Pa*s'),[None,um.degK,um.Pa])
    def visc(self,t,p=None):
        if p==None:
            p=self.__pabs.to('Pa').magnitude
        return self.__vnvisc(t,p,self.__abs.magnitude)
    
    def Pr(self,t,p=None):
        if p==None:
            p=self.__pabs
            
        c=self.cp(t,p)
        v=self.visc(t,p)
        con=self.cond(t,p)
        
        return (c*v/con).to(um.dimensionless)

        
if __name__=='__main__':
    """
    cod de invatare 
    ---------------
    
    def dens(tmp,press):
        return fld.Props('D','T',tmp+273.15,'P',press,'Water')
        
    udens=un.wrap(dens)
    print udens(un.ufloat(10,0.02),un.ufloat(200,0.01))
    vudens=np.vectorize(udens)
    temps=np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)])
    print temps
    press=np.array([un.ufloat(200,0.02),un.ufloat(200,0.02)])    
    print press
    densits=vudens(temps,press)
    print densits,type(densits)
    print vudens(10,200)
    """
    #-------------Water tests----------------------------------
    print 'Running Water tests.......'    
    water=Water()
    #-------test temp is of type: ndarray of UFloat -----------
    temps=np.array([un.ufloat(10,0.02),un.ufloat(20,0.02)])
    dens=water.densit(Q_(temps,um.degC))
    print 'dens=',dens  
    #out:   
    #dens= [999.7496259274476+/-0.0017598342895507812
    #       998.2523477831393+/-0.004131866455078125]
    
    #-------test temp is of type: UFloat-----------------------
    temp=Q_(un.ufloat(10,0.02),um.degC)
    print 'dens=',water.densit(temp)
    #out:
    #dens= 999.7496+/-0.0018    
    
    #------test temp is of type: float-------------------------
    temp=Q_(10.0,um.degC)
    print 'dens=',water.densit(temp)
    #out:
    #dens= 999.749625927
        
    #------test Pr temp is of type ndarray of float------------    
    print 'Pr=',water.Pr(Q_(temps,um.degC))
    #out: (trebuie verificate)
    #Pr= [9.463021047858366+/-0.006206084509439321
    #     7.006353702684516+/-0.0038692636351550565]
    
    #--------------Air tests----------------------------------
    print 'Running Air tests......'
    aer=Air(Rh=0.4)
    
    print 'dens=',aer.densit(Q_(30,um.degC))        
    print 'Pr=',aer.Pr(Q_(40,um.degC))
    print 'Pr=',aer.Pr(Q_(un.ufloat(40,0.1),um.degC))
    print 'Pr=',aer.Pr(Q_([un.ufloat(20,0.1),un.ufloat(40,0.2)],um.degC))
    
    
    
    
        
        
        
        
        