# -*- coding: utf-8 -*-
"""
Created on Wed May 28 08:37:31 2014

@author: Vlad
"""
from __future__ import division
import pint
from __init__ import um,Q_
from excel import OpenXls,InstrErr
import Conversion as conv
from incert import ucreate
import numpy as np
import scipy as sp
import uncertainties as un
import uncertainties.unumpy as unp
import Fluide as flds
import time
#import Structure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ploting import uplot,unitsPrnt
from scipy.optimize import brentq
from unequal import uequal
import dill
from joblib import Parallel, delayed
import multiprocessing
 
def balanceEq(data):
    if not isinstance(data,OpenXls):
        Exception("data trebuie sa fie OpenXls")
    return None


@um.wraps(um.degK,[um.degK,um.degK])
def Delta(x,y):
    return x-y

def Reynolds(v,dh,T,fld):
    if not isinstance(fld,flds.Water) or not isinstance(fld,flds.Air):
        Exception('fld trebuie sa fie ori apa ori aer')
        
    return fld.densit(T)*v*dh/fld.visc(T)

def uflequal(a,b):
    if isinstance(a,float) or isinstance(a,int):
        au=un.ufloat(a,0)
    else:
        au=a
    if isinstance(b,float) or isinstance(b,int):
        bu=un.ufloat(b,0)
    else:
        bu=b
        
    if au.n>bu.n:
        s=bu
        bu=au
        au=s
    aint=[au.n-au.s,au.n+au.s]
    bint=[bu.n-bu.s,au.n+bu.s]
    if (min(aint[1],bint[1])-max(aint[0],bint[0]))>=0:
        return True
    else:
        return False

def indfilter(a,val):
    try:
        x=a.to(val.units)
    except pint.DimensionalityError:
        return None
   
    y=np.ones(a.magnitude.shape)*val.magnitude
    x=x.magnitude
    
def unique(a):
    if isinstance(a,np.ndarray):
        x=[k.magnitude for k in a]
    else:
        x=a.magnitude
    xu=[k.n for k in x]    
    u,ind=np.unique(xu,return_index=True)
    return a[ind]

@um.wraps(um.dimensionless,[um.dimensionless,um.dimensionless],False)   
def Eps(ntu,mu):
    """
    cele doua fluide sunt neamestecate
    """
    return 1.-unp.exp(1/mu*unp.pow(ntu,0.22)*(unp.exp(-mu*unp.pow(ntu,0.78))-1))

@um.wraps(um.dimensionless,[um.dimensionless,um.dimensionless],False)
@np.vectorize
@un.wrap
def NTU(eps,mu):
    """
    inversul functiei Eps
    """
    def f(ntu):
        return eps-(1.-np.exp(1/mu*np.power(ntu,0.22)*(np.exp(-mu*np.power(ntu,0.78))-1)))
    #trebuie vectorizata functia NTU din cauza functiei brentq    
    ret=brentq(f,0,100)
    return ret
    
@um.wraps(um.dimensionless,[um.dimensionless,um.dimensionless],False)
def crit_gni(Re,Pr):
    """
    Calculeaza criteriul Gnielinski
    
    Intrare
    -------
    Re: criteriul Reynolds
    Pr: criteriul Prandtl
    
    Iesire
    ------
    Nu: criteriul Nusselt
    """
    
    fric=(0.79*unp.log(Re)-1.64)**(-2)
    return (fric/8.)*(Re-1000.)*Pr/(1+12.7*unp.sqrt(fric/8.)*(Pr**(2./3)-1))

@um.wraps(um.dimensionless,[um.kg/um.m**3,um.m/um.s,um.milliH2O,um.m**2,um.m**2,um.dimensionless,um.dimensionless],False)
@np.vectorize
@un.wrap
def crit_friction(densit,vit,dp,Ac,At,Kc,Ke):
    return Ac/At*(2*dp/(densit*vit**2)-Ke-Kc)
 
@um.wraps(um.parse_expression('W/(m**2*K)'),[um.W/um.K,um.m**2,um.m**2,um.m,um.m],False)   
@np.vectorize
@un.wrap
def detconvect(kA,Af,A,h,g):
    
    def ka_func(alpha):
        ml=np.sqrt(2*alpha/(237.0*g))*h/2
        eta=1-Af/A*(1-np.tanh(ml)/ml)
        return kA-eta*A*alpha
        
    ret=brentq(ka_func,0.001,1000)
    return ret

@um.wraps(um.parse_expression('W/(m**2*K)'),[um.W/um.K,um.m**2,um.m**2,um.m,um.m],False)
def paralleldetconvect(kA,Af,A,h,g):
    
    def convect(ka,af,a,inalt,gros):

        def ka_func(alpha):
            ml=np.sqrt(2*alpha/(237.0*gros))*inalt/2
            eta=1-af/a*(1-np.tanh(ml)/ml)
            return kA-eta*A*alpha
        
        ret=brentq(ka_func,0.001,1000)
        return ret

    num_core=multiprocessing.cpu_count()
    ret = Parallel(num_core-1)(delayed(convect)(kA[i],Af,A,h,g) for i in range(len(kA)))
    return ret
 
class Radiator():
    def __init__(self,Name,fileName,**kwargs):
        self.update(Name,fileName,**kwargs)
        
    def update(self,Name,fileName,**kwargs):
        """
        se adauga in timp ce ne trebuie
        """
        if kwargs!=None:
            if 'figoffset' in kwargs.keys():
                self.figoffset=kwargs['figoffset']
            else:
                self.figoffset=0
        self.Nume=Name
        #------------------------------------------------------------------
        #            incarcam fisierul cu datele de test
        #------------------------------------------------------------------
        if kwargs==None:
            self.data=OpenXls(fileName)
        else:
            if 'lastrow' in kwargs.keys():
                self.data=OpenXls(fileName,lastrow=kwargs['lastrow'])
            else:
                self.data=OpenXls(fileName)
        
        #------------------------------------------------------------------
        #             calculam conditiile atmosferice
        #------------------------------------------------------------------        
        self.rh_a = self.data['rh_a'].sum()/len(self.data['rh_a'])/100.0
        self.pb_a = (self.data['pb_a'].sum()/len(self.data['pb_a']))
        
        #self.tab_a = (self.data['tin_a'].sum()/len(self.data['tin_a']))
        self.tab_a=Q_(np.array([t.magnitude for t in self.data['tin_a']]).sum()/len(self.data['tin_a']),self.data['tin_a'][0].units)
        #-----------------------------------------------------------------
        #             initializam fluidele
        #-----------------------------------------------------------------
        self.water=flds.Water()
        self.air=flds.Air(self.tab_a,self.pb_a,self.rh_a)
        
    def Calculate(self):
        
        tm_w=((self.data['tin_w'].to('K')+self.data['tout_w'].to('K'))/2.0).to('degC')
        tm_a=((self.data['tin_a'].to('K')+self.data['tout_a'].to('K'))/2.0).to('degC')

        Qw=self.data['debit_w']*self.water.densit(self.data['tin_w'])*self.water.cp(tm_w)*Delta(self.data['tin_w'],self.data['tout_w'])
        Qa=self.data['debit_a']*self.air.cp(tm_a)*Delta(self.data['tout_a'],
                self.data['tin_a'])        
        
        Qm=(Qw+Qa)/2.0

        #-------------------------------------------------------------------
        #   calculam eps=(ties_a-tin_a)/(tin_w-tin_a)
        #            mu=m_a*cp_a/(m_w*cp_w)
        #            NTU =f^-1(eps,mu) -se calculeaza cu inversa functiei
        #-------------------------------------------------------------------
        eps=(self.data['tout_a'].to('K')-self.data['tin_a'].to('K'))/(self.data['tin_w'].to('K')-self.data['tin_a'].to('K'))
        mu=self.data['debit_a']*self.air.cp(tm_a)/(self.data['debit_w']*self.water.densit(self.data['tin_w'])*self.water.cp(tm_w))
        mu=mu.to(um.dimensionless)
        ntu=NTU(eps,mu)
        
        #-------------------------------------------------------------------
        #                       kA_total=NTU*(m_a*cp_a)
        #-------------------------------------------------------------------
        kA=ntu*(self.data['debit_a']*self.air.cp(tm_a))
        kA=kA.to('W/K')
        
        #-------------------------------------------------------------------
        #          calculam vitezele si Re pentru ambele fluide
        #-------------------------------------------------------------------
        vit_w=(self.data['debit_w']/self.data.Ac_w).to('m/s')
        vit_a=self.data['debit_a']/(self.air.densit(tm_a)*self.data.Ac_a)
    
        Re_w=Reynolds(vit_w,self.data.Dh_w,tm_w,self.water)
        Re_w=Re_w.to(um.dimensionless)
        Re_a=Reynolds(vit_a,self.data.Dh_a,tm_a,self.air)
        Re_a=Re_a.to(um.dimensionless)
        Pr_a=self.air.Pr(tm_a)
        
        #-------------------------------------------------------------------
        #                 calculam kA_w folosind criteriul Gnielinski
        #-------------------------------------------------------------------
        alpha_w=crit_gni(Re_w,self.water.Pr(tm_w))*self.water.cond(tm_w)/self.data.Dh_w
        alpha_w=alpha_w.to('W/(m**2*K)')
        kA_w=self.data.At_w*alpha_w
        al_cond=237.0*um.parse_expression('W/m/K')
        
        #calculam aria peretelui despartitor necesar pentru eficienta
        # aria de schimb termic a peretelui nu este dubla
        Ap=self.data.Nr_a*(self.data.Lungime*self.data.Grosime) 
        kA_p=(al_cond*Ap)/self.data.PereteDesp
        kA_p=kA_p.to('W/K')
        kA_a=1/(1/kA-1/kA_p-1/kA_w)     
        Af=self.data.At_a-Ap
        
        #--------------------------------------------------------------------
        #           determinam convectia pe aer folosind metode
        #           numerice de determinare vezi detconvect
        #-------------------------------------------------------------------
        #TODO: de modif detconvect pentru executia in paralel a calculului
        alpha_a=detconvect(kA_a,Af,self.data.At_a,self.data.Inalt_a,
                           self.data.Gros_a)

        #alpha_a=paralleldetconvect(kA_a,Af,self.data.At_a,self.data.Inalt_a,
        #                   self.data.Gros_a)
        
        Nu_a=alpha_a*self.data.Dh_a/self.air.cond(tm_a)
        Nu_a=Nu_a.to(um.dimensionless)
        
        J_a=Nu_a/(Re_a*Pr_a**(1/3))
        
        #-------------------------------------------------------------------
        #          calculam coeficientul de frecare pe aripioara
        #-------------------------------------------------------------------
        
        fr_a=crit_friction(self.air.densit(tm_a),vit_a,self.data['dp_a'],
                           self.data.Ac_a,self.data.At_a,0.0,0.0)
        #-------------------------------------------------------------------
        #             adaugam datele de return la clasa 
        #-------------------------------------------------------------------
        self.__dict__.update({'Qm':Qm.to('kW'),
                              'deb_a':self.data['debit_a'],
                              'deb_w':self.data['debit_w'],
                              'tin_a':self.data['tin_a'],
                              'tout_a':self.data['tout_a'],
                              'tin_w':self.data['tin_w'],
                              'tout_w':self.data['tout_w'],
                              'h':self.data.Inalt_a,
                              'p':self.data.Pas_a,
                              'Dh':self.data.Dh_a,
                              'L':self.data.Grosime,
                              'Ac_a':self.data.Ac_a,
                              'Ac_w':self.data.Ac_w,
                              'At_a':self.data.At_a,
                              'At_w':self.data.At_w,
                              'Ap':Ap,
                              'v':vit_a,
                              'v_w':vit_w,
                              'Re':Re_a,
                              'Re_w':Re_w,
                              'Pr':self.air.Pr(tm_a),
                              'Nu':Nu_a,
                              'J':J_a,
                              'fr':fr_a,
                              'alpha':alpha_a,
                              'alpha_w':alpha_w,
                              'kA':kA
                              })   
        self._reduce()
        
        
    def _reduce(self):
        def unique(lst):
            ret=[lst[0]]
            for x in lst:
                ins=False
                for y in ret:
                    if uequal(y,x):
                        ins=True
                        break
                if not ins:
                    ret.append(x)
            return np.array(ret)
       
        Reu=unique(self.Re.magnitude)
        rRe=[]
        rNu=[]
        rfr=[] 
        rPr=[]
        
        for re in Reu:
            inds=[np.where(self.Re.magnitude==x)[0][0] 
                for x in self.Re.magnitude if uequal(re,x)]
            res=self.Re.magnitude[inds]            
            rRe.append(np.mean(res))
            rPr.append(np.mean(self.Pr.magnitude[inds]))
            rNu.append(np.mean(self.Nu.magnitude[inds]))
            rfr.append(np.mean(self.fr.magnitude[inds]))
        
        self.__dict__.update({'rRe':Q_(rRe,um.dimensionless),
                              'rNu':Q_(rNu,um.dimensionless),
                              'rfr':Q_(rfr,um.dimensionless),
                              'rPr':Q_(rPr,um.dimensionless)})
            
    def save_csv(self,*datasave):
        import sys
        import csv
        if not datasave:
            datasave=('Qm','deb_a','tin_a','tout_a','v','Re','alpha','deb_w','tin_w','tout_w','v_w','Re_w','alpha_w')
        csvdta=zip(*(self.__dict__[k] for k in datasave)) 

        with open("%s.csv"%self.Nume,"wb") as file:
            writer=csv.writer(file,delimiter=',')
            writer.writerow(datasave)
            header=tuple(a.units for a in csvdta[0])
            writer.writerow(header)
            for line in csvdta:
                l=tuple(a.magnitude for a in line)
                writer.writerow(l)

            writer.writerow(("At_a",self.At_a.magnitude,self.At_a.units))
            writer.writerow(("Ac_a",self.Ac_a.magnitude,self.Ac_a.units))
            writer.writerow(("At_w",self.At_w.magnitude,self.At_w.units))
            writer.writerow(("Ac_a",self.Ac_w.magnitude,self.Ac_w.units))
    
    def plot_Qm3D(self,save=False):
        fig=plt.figure(self.figoffset)
        ax=Axes3D(fig)
        ax.scatter(strip(self.deb_w),strip(self.deb_a),strip(self.Qm))
        ax.set_xlabel('debit apa [{:}]'.format(unitsPrnt(self.deb_w[0])))
        ax.set_ylabel('debit aer [{:}]'.format(unitsPrnt(self.deb_a[0])))
        ax.set_zlabel('Flux termic [{:}]'.format(unitsPrnt(self.Qm[0])))
        plt.title(u'Performantele {:}'.format(self.Nume))
        if save==True: 
            plt.savefig('{:}_3d.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),
                        orientation='landscape',
                        facecolor='w',
                        edgecolor='k')
                        
        self.figoffset=self.figoffset+1
        
    def plot_Qm(self,save=False):
        wdebite=unique(self.data['debit_w'])
        qm=np.array([self.data['debit_w'],self.data['debit_a'],
                     (self.Qm.to('kW'))])
        qm=qm.T
        
        #plt.figure(self.figoffset+1)
        for i,val in zip(range(len(wdebite)),wdebite):
            #---------filtram datele
            plt.figure(self.figoffset+i)            
            #fig,axes=plt.subplots(2,2)            
            qm1=np.array([s for s in qm if uflequal(s[0],val)])
            uplot(qm1[:,1],qm1[:,2],func=lambda x,a,b:a*x**b,fmt='k.',
                  noerr=True)
            plt.title('{:s}= {:~P}'.format(r'$\dot{m}_w$',val))  
            plt.xlabel('{:s} [{:}]'.format(r'$\dot{m}_a$',unitsPrnt(qm1[0,1])))
            plt.ylabel('{:s} [{:}]'.format(r'$\dot{Q}_a$',unitsPrnt(qm1[0,2])))
            plt.grid(which='both',axis='both')        
            plt.legend(loc='lower right')        
            plt.tight_layout()
            if save==True:        
                plt.savefig('{:}_qm_{:}.png'.format(self.Nume,i),dpi=300,
                            figsize=(166.54/2.54,81/2.54),
                            orientation='landscape',
                            facecolor='w',edgecolor='k')
                        
        self.figoffset=self.figoffset+i+1
    
    def plot_Nu(self,save=False):
        plt.figure(self.figoffset)
        pl,coef=uplot(self.Re,self.Nu,sg=1,func=lambda x,a,b:a*x**b,
                      noerr=True)
        #plt.title(r'%s: $Nu=%.3f\pm %.3f*Re^{%.3f\pm %.3f}$'
        #    %(self.Nume,coef[0].n,coef[0].s,coef[1].n,coef[1].s))
        plt.xlabel('Re')
        plt.ylabel('Nu')
        plt.grid(which='both',axis='both')
        #plt.legend(loc='lower right')
        ax=plt.gca()   
        txt=r'$Nu=%.3f\pm %.3f\cdot Re^{%.3f\pm %.3f}$'%(
            coef[0].n,coef[0].s,coef[1].n,coef[1].s)
        plt.text(0.4,0.05,txt,
                 zorder=5,size=20,bbox=dict(facecolor='aqua', alpha=0.8),
                 transform=ax.transAxes)        
        plt.tight_layout()
        if save==True:
            plt.savefig('{:}_Nu.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),orientation='landscape',
                        facecolor='w',edgecolor='k')
        
        self.figoffset=self.figoffset+1
   
    def plot_J(self,save=False):
        plt.figure(self.figoffset)
        pl,coef=uplot(self.Re,self.J,sg=1,func=lambda x,a,b:a*x**b,noerr=True)
        #plt.title(r'%s: $J=%.3f\pm %.3f*Re^{%.3f\pm %.3f}$'
        #    %(self.Nume,coef[0].n,coef[0].s,coef[1].n,coef[1].s))
        plt.xlabel('Re')
        plt.ylabel('J')
        plt.grid(which='both',axis='both')
        #plt.legend(loc='lower right')
        ax=plt.gca()    
        txt=r'$J=%.3f\pm %.3f\cdot Re^{%.3f\pm %.3f}$'%(coef[0].n,coef[0].s,
                                                    coef[1].n,coef[1].s)
        plt.text(0.4,0.05,txt,
                 zorder=5,size=20,bbox=dict(facecolor='aqua', alpha=0.8),
                 transform=ax.transAxes)
        plt.tight_layout()
        if save==True:
            plt.savefig('{:}_J.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),orientation='landscape',
                        facecolor='w',edgecolor='k')
        
        self.figoffset=self.figoffset+1
        return self.figoffset-1
    
    def plot_f(self,save=False,lim=2600):
        plt.figure(self.figoffset) 
        
        inds_low=[i for i,j in enumerate(self.Re) if j<=lim]
        inds_high=[i for i,j in enumerate(self.Re) if j>=lim]
        
        pl,coef_low=uplot(self.Re[inds_low],self.fr[inds_low],sg=1,
                      func=lambda x,a,b:a*x**b,fmt='k.',
                      noerr=True,slim=lim)
       
        pl,coef_high=uplot(self.Re[inds_high],self.fr[inds_high],sg=1,
                      func=lambda x,a,b:a*x**b,fmt='k.',
                      noerr=True,ilim=lim)
        ax=plt.gca()
        text=r'$c_f=\left \{ \stackrel{%.3f\pm %.3f\cdot Re^{%.3f\pm %.3f}\;Re\leq 2600}{%.3f\pm %.3f\cdot Re^{%.3f\pm %.3f}\;Re>2600}\right. $'%(\
                                coef_low[0].n,coef_low[0].s,\
                                coef_low[1].n,coef_low[1].s,\
                                coef_high[0].n,coef_high[0].s,\
                                coef_high[1].n,coef_high[1].s)
        
        
        #plt.title(r'%s: $c_f=%.3f\pm %.3f*Re^{%.3f\pm %.3f}$'
        #    %(self.Nume,coef[0].n,coef[0].s,coef[1].n,coef[1].s))
        
        plt.xlabel('Re')
        plt.ylabel(r'$c_f$')
        plt.grid(which='both',axis='both')
        
        #plt.legend(loc='upper right')
        
        plt.text(0.3,0.85,text,zorder=5,size=20,bbox=dict(facecolor='aqua', alpha=0.8),
                 transform=ax.transAxes)
        plt.tight_layout()
        if save==True:
            plt.savefig('{:}_cf.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),orientation='landscape',
                        facecolor='w',edgecolor='k')
        
        self.figoffset=self.figoffset+1
        return self.figoffset-1
        
    def plot_rNu(self,save=False,fig=None):
        if fig==None:
            plt.figure(self.figoffset)
        else:
            plt.figure(fig)
            
        pl,coef=uplot(self.rRe,self.rNu,sg=1,func=lambda x,a,b: a*x**b,
                      noerr=True)
        if fig==None:
            plt.title(r'%s: $Nu=%.3f\pm %.3f*Re^{%.3f\pm %.3f}$'
                %(self.Nume,coef[0].n,coef[0].s,coef[1].n,coef[0].s))
            plt.xlabel('Re')
            plt.ylabel('Nu')
            plt.grid(which='both',axis='both')
            #plt.legend(loc='lower right')
            plt.tight_layout()
        
        if fig==None and save:
            plt.savefig('{:}_rNu.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),orientation='landscape',
                        facecolor='w',edgecolor='k')
        
        if fig==None:
            self.figoffset=self.figoffset+1
            return self.figoffset-1
        else:
            return fig
            
    def plot_rf(self,save=False,fig=None):
        if fig==None:
            plt.figure(self.figoffset)
        else:
            plt.figure(fig)
            
        pl,coef=uplot(self.rRe,self.rfr,sg=1,func=lambda x,a,b: a*x**b,
                      noerr=True)
        if fig==None:
            plt.title(r'%s: $C_f=%.3f\pm %.3f*Re^{%.3f\pm %.3f}$'
                %(self.Nume,coef[0].n,coef[0].s,coef[1].n,coef[0].s))
            plt.xlabel('Re')
            plt.ylabel('$C_f$')
            plt.grid(which='both',axis='both')
            plt.legend(loc='lower right')
            plt.tight_layout()
        
        if fig==None and save:
            plt.savefig('{:}_rNu.png'.format(self.Nume),dpi=300,
                        figsize=(166.54/2.54,81/2.54),orientation='landscape',
                        facecolor='w',edgecolor='k')
        
        if fig==None:
            self.figoffset=self.figoffset+1
            return self.figoffset-1
        else:
            return fig
            

        
def strip(ar):
    return unp.nominal_values(ar.magnitude)


        
def testLoadSaveRa():
    ra=Radiator('Ra28788-0',
                'D:/Documents/Lucrare Doctorat/Experimental/RA 28788-0.xls',
                figoffset=0)
    ra.Calculate()
    print 'saving .....'
    with open('{:s}.pkl'.format(ra.Nume),'wb') as output:
        dill.dump(ra,output)
        
        
    print 'loading.....'
    raload=None
    with open('{:s}.pkl'.format(ra.Nume),'rb') as intput:
        raload=dill.load(intput)
        
    print ra.Re[1:3]
    print raload.Re[1:3]
    
def testRa():
    import csv
    import sys
    ra=Radiator('Ra28788-0',
                'E:/Documents/Lucrare Doctorat/Experimental/RA 28788-0.xls',
                figoffset=0)
    print 'calculeaza ......'
    ra.Calculate()
    
    print 'salveaza csv....'
    ra.save_csv()

    print 'plot.......'
    ra.plot_Qm3D()
    print ra.figoffset
    ra.plot_Qm()
    print ra.figoffset
    ra.plot_Nu()
    print ra.figoffset
    ra.plot_J()
    print ra.figoffset
    ra.plot_f(save=True)
    print ra.figoffset
    
    ra.plot_rNu()   
    ra.plot_rf()
    
    

def testbalance():
    start=time.time()
    file_location='D:/Documents/Lucrare Doctorat/Experimental/RA 28788-0.xls'
    data=OpenXls(file_location)
    print 'open time',time.time()-start
    start=time.time()
    rh_a = data['rh_a'].sum()/len(data['rh_a'])/100.0
    pb_a = (data['pb_a'].sum()/len(data['pb_a']))
    tab_a = (data['tin_a'].sum()/len(data['tin_a']))
    end=time.time()
    print 'reading time',end-start
    start=time.time()
    water=flds.Water()
    air=flds.Air(tab_a,pb_a,rh_a)
    tm_w=(data['tin_w']+data['tout_w'])/2.0
    tm_a=(data['tin_a']+data['tout_a'])/2.0
    Qw=data['debit_w']*water.densit(data['tin_w'])*water.cp(tm_w)*Delta(data['tin_w'],data['tout_w'])
    Qa=data['debit_a']*air.cp(tm_a)*Delta(data['tout_a'],data['tin_a'])
    
    Qm=(Qw+Qa)/2.0

    eps=(data['tout_a'].to('K')-data['tin_a'].to('K'))/(data['tin_w'].to('K')-data['tin_a'].to('K'))
    mu=data['debit_a']*air.cp(tm_a)/(data['debit_w']*water.densit(data['tin_w'])*water.cp(tm_w))
    mu=mu.to(um.dimensionless)
    ntu=NTU(eps,mu)
    kA=ntu*(data['debit_a']*air.cp(tm_a))
    kA=kA.to('W/K')
    
    vit_w=data['debit_w']/data.Ac_w
    vit_a=data['debit_a']/(air.densit(tm_a)*data.Ac_a)
    
    Re_w=Reynolds(vit_w,data.Dh_w,tm_w,water)
    Re_w=Re_w.to(um.dimensionless)
    Re_a=Reynolds(vit_a,data.Dh_a,tm_a,air)
    Re_a=Re_a.to(um.dimensionless)
    
    #calculam kA_w
    alpha_w=crit_gni(Re_w,water.Pr(tm_w))*water.cond(tm_w)/data.Dh_w
    alpha_w=alpha_w.to('W/(m**2*K)')
    kA_w=data.At_w*alpha_w
    al_cond=237.0*um.parse_expression('W/m/K')
    Ap=2*data.Nr_a*(data.Lungime*data.Grosime)
    kA_p=(al_cond*Ap)/data.PereteDesp
    kA_p=kA_p.to('W/K')
    kA_a=1/(1/kA-1/kA_p-1/kA_w)     
    Af=data.At_a-Ap
    alpha_a=detconvect(kA_a,Af,data.At_a,data.Inalt_a,data.Gros_a)
    Nu_a=alpha_a*data.Dh_a/air.cond(tm_a)
    Nu_a=Nu_a.to(um.dimensionless)
        
    end=time.time()
    print 'calculation time',end-start
    #-----construim un tabel cu debitele necesare pentru plotare
    qm=np.array([data['debit_w'],data['debit_a'],(Qm.to('kW'))])
    qm=qm.T
    #ploting(qm,data)
    plt.figure(10)
    pl,coefs=uplot(Re_a,Nu_a,func=lambda re,a,b:a*re**b)
    plt.tight_layout()
    plt.show()
    for coef in coefs:
        print '{:P}'.format(coef)

def ploting(qm,data):
    print 'ploting.....'
    wdebite=unique(data['debit_w'])
    for i,val in zip(range(len(wdebite)),wdebite):
        #---------filtram datele
        plt.figure(i)        
        qm1=np.array([s for s in qm if uflequal(s[0],val)])
        uplot(qm1[:,1],qm1[:,2],func=lambda x,a,b:a*x**b,fmt='k.')
        plt.title('Debitul de apa= {:~P}'.format(val))  
        plt.xlabel('debit aer [{:}]'.format(unitsPrnt(qm1[0,1])))
        plt.ylabel('flux mediu [{:}]'.format(unitsPrnt(qm1[0,2])))
        plt.grid(which='both',axis='both')        
        plt.legend(loc='lower right')        
        plt.tight_layout()
        plt.savefig('qm_{:}'.format(i),dpi=300,figsize=(166.54/2.54,81/2.54),orientation='landscape',facecolor='w',edgecolor='k')
        plt.show()        
        if i==4:
            print qm1[3:5,2]
    
    
if __name__=='__main__':
    
    
    plt.close('all')
    #testbalance()    
    #struct=Struct(2.543e-3,8.36e-3,5.465,3.25e-3)
    #print struct
    #print struct.Ac
    testRa()
    #testLoadSaveRa()
    