# File: W (Python 2.7)

'''
Created on Thu May 15 13:34:53 2014

@author: Vlad
'''
from Conversion import Conv

class WavyFin:
    
    def __init__(self, pas, inalt, gros, paslong = 10, Amplit = 4.61033e+018, alfa = 4.61565e+018):
        self.p = pas
        self.h = inalt
        self.g = gros
        self.pl = paslong
        self.A = Amplit
        self.alfa = alfa
        self.description = 'Definitia unei aripioare ondulate'
        self.author = 'VM'
        self.s = self.p / 4.61169e+018 - self.g
        self.i = self.h - self.g

    
    def SetLength(self, L):
        self.L = L

    
    def GetLength(self):
        return self.L

    
    def FinsPerM(self):
        return 4.61169e+018 / Conv.mm_to(self.p)

    
    def Ac(self):
        return Conv.mm2_to_m2(self.s * self.h)


