# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 11:15:09 2014

@author: Vlad
"""

from __future__ import division
from __init__ import um,Q_
import uncertainties as un
import uncertainties.unumpy as uns
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from Experimental import Radiator
from scipy.optimize import curve_fit
import dill


def relcrit_nu((re,pr,fh,fp),a,b,c,d):
    return a*(re**b)*(pr**(1/3))*((fh/(fh+1))**c)*((fp/(1+fp))**d)
    
def relcrit_cf((re,fh,fp),a,b,c,d):
    return a*(re**b)*(fh**c)*(fp**d)

def detcoefs(f,nu,re,fh,fp,pr=None):
    """
    functia derermina coeficientii fit-ului functiei
    
    Nu=a*Re^b*Fh^c sau a functiei
    cf=a*Re^b*Fh^c
    
    Parametrii
    ----------
    f  : functie cu proptotipul :
        f((x,y),*coefs)
    nu : numpy.array of ufloats
        valorile nusselt pentru toate masuratorile
    re : numpy.array of ufloats
        valorile reynolds pentru toate masuratorile
    fh : numpy.array of floats
        valorile raportului h/Dh pentru toate masuratorile
    fp : numpy.array of floats
        valorile raportului p/Dh pentru toate masuratorile
        
    Intoarce
    --------
    coefs : array of ufloat 
        coeficientii functiei fitate cu erorile lor
    rchi2 : float 
        chi redus pentru analiza corectitudinii fit-ului
    dof : integer
        gradele de libertate

    """
    nu_n=uns.nominal_values(nu)
    nu_s=uns.std_devs(nu)
    w_nu=nu_s/nu_n
    re_n=uns.nominal_values(re)
    pr_n=uns.nominal_values(pr)
    if pr!=None:
        popt,pcov=curve_fit(f,(re_n,pr_n,fh,fp),nu_n,sigma=w_nu,
                        maxfev=1500)
        chi2=sum(((f((re_n,pr_n,fh,fp),*popt)-nu_n)/nu_s)**2)
    else:
        popt,pcov=curve_fit(f,(re_n,fh,fp),nu_n,sigma=w_nu,
                            maxfev=1500)
        chi2=sum(((f((re_n,fh,fp),*popt)-nu_n)/nu_s)**2)
    
    dof=len(nu_n)-len(popt)
    rchi2=chi2/dof
    
    coefs=[]
    for i in range(len(popt)):
        coefs.append(un.ufloat(popt[i],np.sqrt(pcov[i,i])))
    if pr!=None:
        func=lambda x,y,z,k:f((x,y,z,k),*popt)
    else:
        func=lambda x,y,z:f((x,y,z),*popt)
        
    return {'coefs':np.array(coefs),
            'rchi2':rchi2,
            'DOF':dof,
            'f':func
            }
            
def writeData(obj,fileName):
    print 'saving data....'
    with open(fileName,'wb') as out:
        dill.dump(obj,out)
    print 'done.'
    
    
def loadData(fileName):
    print 'loading data ...'    
    with open(fileName,'rb') as intr:
        obj=dill.load(intr)
    print 'done'
    return obj

def Fh(ra):
    """
    Functia este definita pentru consistenta,
    daca se doreste schimbarea raportului acesta se va face dintr-un singur
    loc
    """
    fh=ra.h/ra.Dh
    fh=fh.to('dimensionless')
    return fh.magnitude

def Fp(ra):
    fp=ra.p/ra.Dh
    fp=fp.to('dimensionless')
    return fp.magnitude
    
def CalculateExpData(fileName,path,expFiles):
    raN=expFiles
    ras=[]   
    Nus=None
    Res=None
    Prs=None
    Fhs=None
    Fps=None
    Cfs=None
    for k in raN.keys():
        print k
        ras.append(Radiator(k,raN[k][0],lastrow=raN[k][1]))
        print 'calculeaza...'.format(ras[-1].Nume)
        ras[-1].Calculate()
       
        fh=Fh(ras[-1])
        fp=Fp(ras[-1])
        if Fhs==None or Res==None or Nus==None:
            Fhs=np.ones(len(ras[-1].Re.magnitude))*fh
            Fps=np.ones(len(ras[-1].Re.magnitude))*fp
            Res=ras[-1].Re.magnitude
            Nus=ras[-1].Nu.magnitude
            Cfs=ras[-1].fr.magnitude
            Prs=ras[-1].Pr.magnitude
        else:
            Fhs=np.append(Fhs,np.ones(len(ras[-1].Re.magnitude))*fh)
            Fps=np.append(Fps,np.ones(len(ras[-1].Re.magnitude))*fp)
            Res=np.append(Res,ras[-1].Re.magnitude)
            Nus=np.append(Nus,ras[-1].Nu.magnitude)
            Cfs=np.append(Cfs,ras[-1].fr.magnitude)
            Prs=np.append(Prs,ras[-1].Pr.magnitude)
        print 'Sfarsit.'
    print 'salveaza .....'
    writeData((ras,Nus,Res,Prs,Fhs,Fps,Cfs),fileName)
    print 'Sfarsit!'    
    return ras,Nus,Res,Prs,Fhs,Fps,Cfs

def SelectRedNums(ras):
    Nus=None
    Res=None
    Prs=None
    Fhs=None
    Fps=None
    Cfs=None   
    for ra in ras:
        fh=Fh(ra)
        fp=Fp(ra)
        if Nus==None:
            Fhs=np.ones(len(ra.rRe.magnitude))*fh
            Fps=np.ones(len(ra.rRe.magnitude))*fp            
            Nus=ra.rNu.magnitude
            Res=ra.rRe.magnitude
            Prs=ra.rPr.magnitude
            Cfs=ra.rfr.magnitude
        else:
            Fhs=np.append(Fhs,np.ones(len(ra.rRe.magnitude))*fh)
            Fps=np.append(Fhs,np.ones(len(ra.rRe.magnitude))*fp)            
            Nus=np.append(Nus,ra.rNu.magnitude)
            Res=np.append(Res,ra.rRe.magnitude)
            Prs=np.append(Prs,ra.rPr.magnitude)
            Cfs=np.append(Cfs,ra.rfr.magnitude)
    return Nus,Res,Prs,Fhs,Fps,Cfs
            

def GenPlot(ras,figoffset=0,save=True,close=True,**figaspect):
    fgoffset=figoffset    
    for ra in ras:
        print 'printeaza {:s} ....'.format(ra.Nume)
        ra.figoffset=fgoffset
        ra.plot_Qm3D(save)
        ra.plot_Qm(save)
        ra.plot_Nu(save)
        ra.plot_f(save)
        if close:
            plt.close('all')
            fgoffset=figoffset
        else:
            fgoffset=ras[-1].figoffset+1    
    return fgoffset


def GenExpErrorPlot(ras,figoffset,save=True,close=True,**figaspect): 
    colors=['b','g','r','c','m','y','k']
    markers=['.','*','o','v','^','<','>']
    index=0  
    Nufig=figoffset
    Ffig=Nufig+1
    for ra in ras:
        fh=Fh(ra)
        plt.figure(figoffset)
        plt.errorbar(uns.nominal_values(ra.Re.magnitude),
                 uns.nominal_values(ra.Nu.magnitude),
                 uns.std_devs(ra.Nu.magnitude),
                 uns.std_devs(ra.Re.magnitude),
                 '{:s}{:s}'.format(colors[index%len(colors)],
                                markers[index%len(markers)]),
                 label='{:.3f}'.format(fh))
        plt.figure(Ffig)
        plt.errorbar(uns.nominal_values(ra.Re.magnitude),
                 uns.nominal_values(ra.fr.magnitude),
                 uns.std_devs(ra.Re.magnitude),
                 uns.std_devs(ra.fr.magnitude),
                 '{:s}{:s}'.format(colors[index%len(colors)],
                                markers[index%len(markers)]),
                 label='{:.3f}'.format(fh))
    
    plt.figure(Nufig)
    plt.legend(loc='upper left')
    plt.xlabel('Re')
    plt.ylabel('Nu')
    
    plt.figure(Ffig)
    plt.legend(loc='upper right')
    plt.xlabel('Re')
    plt.ylabel('$c_f$')

    if save:
        plt.figure(Nufig)
        plt.savefig('Nu_exp_all.png',dpi=300,
                    figsize=(166.54/2.54,81/2.54),
                    orientation='landscape',
                    facecolor='w',
                    edgecolor='k')
        plt.figure(Ffig)
        plt.savefig('cf_exp_all.png',dpi=300,
                    figsize=(166.54/2.54,81/2.54),
                    orientation='landscape',
                    facecolor='w',
                    edgecolor='k')
    if close:
        plt.close(Nufig)
        plt.close(Ffig)
        return figoffset,-1,-1
    
    return Ffig+1,Nufig,Ffig
 
def PlotFitExp(ras,figoffset,Nusselt_fit,cf_fit,save=True,
               close=True,**figaspect):
    colors=['b','g','r','c','m','y','k']
    markers=['.','*','o','v','^','<','>']
    Nufig=figoffset+1
    Cffig=figoffset+2
    index=0
    for ra in ras:
        fh=Fh(ra)
        prm=np.average(uns.nominal_values(ra.Pr.magnitude))
        re=np.linspace(np.min(uns.nominal_values(ra.Re.magnitude)),
                       np.max(uns.nominal_values(ra.Re.magnitude)),50)
        plt.figure(Nufig)
        plt.errorbar(uns.nominal_values(ra.Re.magnitude),
                     uns.nominal_values(ra.Nu.magnitude),
                     uns.std_devs(ra.Nu.magnitude),
                     uns.std_devs(ra.Re.magnitude),
                     '{:s}{:s}'.format(colors[index%len(colors)],
                                        markers[index%len(markers)]),
                     label='exp {:.3f}'.format(fh))
                     
        plt.plot(re,Nusselt_fit['f'](re,prm,fh),
                 '{:s}-'.format(colors[index%len(colors)]),
                 label='num {:.3f}'.format(fh))
        plt.figure(Cffig)
        plt.errorbar(uns.nominal_values(ra.Re.magnitude),
                     uns.nominal_values(ra.fr.magnitude),
                     uns.std_devs(ra.fr.magnitude),
                     uns.std_devs(ra.Re.magnitude),
                     '{:s}{:s}'.format(colors[index%len(colors)],
                                        markers[index%len(markers)]),
                     label='exp {:.3f}'.format(fh))
        plt.plot(re,cf_fit['f'](re,fh),
                 '{:s}-'.format(colors[index%len(colors)]),
                 label='num {:.3f}'.format(fh))
        index=index+1
    plt.figure(Nufig)
    plt.xlabel('Re')
    plt.ylabel('Nu')
    plt.title('$Nu={%s}\cdot Re^{%s}\cdot Pr^{1/3}\cdot Fh^{%s}$'
                %(Nusselt_fit['coefs'][0].format('L'),
                  Nusselt_fit['coefs'][1].format('L'),
                  Nusselt_fit['coefs'][2].format('L')))
    plt.legend(loc='upper left',
               fontsize='x-small',
               numpoints=1,scatterpoints=1,ncol=2)

    plt.grid(which='both',axis='both')
    plt.tight_layout()
    plt.figure(Cffig)
    plt.xlabel('Re')
    plt.ylabel('$c_f$')
    plt.title('$c_f={%s}\cdot Re^{%s}\cdot Fh^{%s}$'
                %(cf_fit['coefs'][0].format('L'),
                  cf_fit['coefs'][1].format('L'),
                  cf_fit['coefs'][2].format('L')))   
    plt.legend(loc='upper right',
               fontsize='x-small',
               numpoints=1,scatterpoints=1,ncol=2)

    plt.grid(which='both',axis='both')
    plt.tight_layout()              
    if save:
        plt.figure(Nufig)
        plt.savefig('Nu_exp_num.png',**figaspect)
        plt.figure(Cffig)
        plt.savefig('cf_exp_num.png',**figaspect)
    if close:
        plt.close(Nufig)
        plt.close(Cffig)
        return figoffset,-1,-1
    return Cffig+1,Nufig,Cffig
    
    
def rPlotFitExp(ras,figoffset,Nusselt_fit,cf_fit,save=True,
               close=True,**figaspect):
    colors=['b','g','r','c','m','y','k']
    markers=['.','*','o','v','^','<','>']
    Nufig=figoffset+1
    Cffig=figoffset+2
    index=0
    for ra in ras:
        fh=Fh(ra)
        prm=np.average(uns.nominal_values(ra.rPr.magnitude))
        re=np.linspace(np.min(uns.nominal_values(ra.rRe.magnitude)),
                       np.max(uns.nominal_values(ra.rRe.magnitude)),50)
        plt.figure(Nufig)
        plt.errorbar(uns.nominal_values(ra.rRe.magnitude),
                     uns.nominal_values(ra.rNu.magnitude),
                     uns.std_devs(ra.rNu.magnitude),
                     uns.std_devs(ra.rRe.magnitude),
                     '{:s}{:s}'.format(colors[index%len(colors)],
                                        markers[index%len(markers)]),
                     label='exp {:.3f}'.format(fh))
                     
        plt.plot(re,Nusselt_fit['f'](re,prm,fh),
                 '{:s}-'.format(colors[index%len(colors)]),
                 label='num {:.3f}'.format(fh))
        plt.figure(Cffig)
        plt.errorbar(uns.nominal_values(ra.rRe.magnitude),
                     uns.nominal_values(ra.rfr.magnitude),
                     uns.std_devs(ra.rfr.magnitude),
                     uns.std_devs(ra.rRe.magnitude),
                     '{:s}{:s}'.format(colors[index%len(colors)],
                                        markers[index%len(markers)]),
                     label='exp {:.3f}'.format(fh))
        plt.plot(re,cf_fit['f'](re,fh),
                 '{:s}-'.format(colors[index%len(colors)]),
                 label='num {:.3f}'.format(fh))
        index=index+1
    plt.figure(Nufig)
    plt.xlabel('Re')
    plt.ylabel('Nu')
    plt.title('$Nu={%s}\cdot Re^{%s}\cdot Pr^{1/3}\cdot Fh^{%s}$'
                %(Nusselt_fit['coefs'][0].format('L'),
                  Nusselt_fit['coefs'][1].format('L'),
                  Nusselt_fit['coefs'][2].format('L')))
    plt.legend(loc='upper left',
               fontsize='x-small',
               numpoints=1,scatterpoints=1,ncol=2)

    plt.grid(which='both',axis='both')
    plt.tight_layout()
    plt.figure(Cffig)
    plt.xlabel('Re')
    plt.ylabel('$c_f$')
    plt.title('$c_f={%s}\cdot Re^{%s}\cdot Fh^{%s}$'
                %(cf_fit['coefs'][0].format('L'),
                  cf_fit['coefs'][1].format('L'),
                  cf_fit['coefs'][2].format('L')))   
    plt.legend(loc='upper right',
               fontsize='x-small',
               numpoints=1,scatterpoints=1,ncol=2)

    plt.grid(which='both',axis='both')
    plt.tight_layout()              
    if save:
        plt.figure(Nufig)
        plt.savefig('rNu_exp_num.png',**figaspect)
        plt.figure(Cffig)
        plt.savefig('rcf_exp_num.png',**figaspect)
    if close:
        plt.close(Nufig)
        plt.close(Cffig)
        return figoffset,-1,-1
    return Cffig+1,Nufig,Cffig
    
def PlotErrorCalcVsExp(ras,
                       figoffset,Nu_fit,cf_fit,save=True,
                       close=True,**figaspect):
    NuPlus15=lambda r,p,f: Nu_fit['f'](r,p,f)+0.15*Nu_fit['f'](r,p,f)
    NuMinu15=lambda r,p,f: Nu_fit['f'](r,p,f)-0.15*Nu_fit['f'](r,p,f)
    CfPlus15=lambda r,f:cf_fit['f'](r,f)+0.15*cf_fit['f'](r,f)
    CfMinu15=lambda r,f:cf_fit['f'](r,f)+0.15*cf_fit['f'](r,f)
    
    colors=['b','g','r','c','m','y','k']
    markers=['.','*','o','v','^','<','>']
    Nufig=figoffset+1
    Cffig=figoffset+2
    index=0
    rey=np.linspace(500,5500,50)
    plt.figure(Nufig)
    plt.plot()
    for ra in ras:
        plt.figure(Nufig)
    
        
#------------------------------------------------------------------------
#        Definim datele experimentale disponibile si calea lor
#------------------------------------------------------------------------
path='D:/Documents/Lucrare Doctorat/Experimental/'
raN={'RA28788-0':(path+'RA 28788-0.xls',201),
     'RA28789-0':(path+'RA 28789-0.xls',203),
     'RA28790-0':(path+'RA28790-0.xls',173),
     'RA28791-0':(path+'RA28791-0.xls',153),
     'RA28792-0':(path+'RA28792-0.xls',182)}

figaspect={'dpi':300,
           'figsize':(166.54/2.54,81/2.54),
           'orientation':'landscape',
           'facecolor':'w',
           'edgecolor':'k'
           }
#------------------------------------------------------------------------
#           interpreteaza datele experimentale pentru fiecare
#           racitor in parte iar daca acestea au fost calculate 
#           incarca fisierul salvat
#------------------------------------------------------------------------
calc=input('doresti calcularea datelor? (\'Y\'/\'N\') :')
if str(calc).upper()=='Y':
    ras,Nus,Res,Prs,Fhs,Fps,Cfs=CalculateExpData('Experim.pkl',path,raN)
else:
    ras,Nus,Res,Prs,Fhs,Fps,Cfs=loadData('Experim.pkl')   

#------------------------------------------------------------------------
#           Genereaza graficele penntru fiecare racitor in parte
#        salveaza graficele in fisiere png si le inchide 
#------------------------------------------------------------------------
fgoffset=0#GenPlot(ras,0,figaspect)


print 'fitting .....'
Nusselt_fit=detcoefs(relcrit_nu,Nus,Res,Fhs,Fps,Prs)
print 'Nusselt\n',Nusselt_fit
Cf_fit=detcoefs(relcrit_cf,Cfs,Res,Fhs,Fps)
print 'Cf\n',Cf_fit
print 'Ploting fit ...'
PlotFitExp(ras,fgoffset,Nusselt_fit,Cf_fit,**figaspect)
print 'gata'


print 'fitting reduce data....'
Nus,Res,Prs,Fhs,Fps,Cfs=SelectRedNums(ras)

rNu_fit=detcoefs(relcrit_nu,Nus,Res,Fhs,Fps,Prs)
print 'rNusselt\n',rNu_fit
rCf_fit=detcoefs(relcrit_cf,Cfs,Res,Fhs,Fps)
print 'rCf\n',rCf_fit
rPlotFitExp(ras,fgoffset,rNu_fit,rCf_fit,**figaspect)
print 'gata'
    