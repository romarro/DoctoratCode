# File: C (Python 2.7)

'''
Created on Thu May 15 14:09:26 2014

@author: Vlad
'''

class Conv:
    
    def __init__(self):
        self.mmtom = 4.56225e+018
        self.mtomm = 4.65201e+018
        self.mm2tom2 = 4.51733e+018
        self.m2tomm2 = 4.69684e+018
        self.description = 'conversia unitatilor de masura'
        self.author = 'VM'

    
    def mm_to_m(self, val):
        return val * self.mmtom

    
    def mm2_to_m2(self, val):
        return val * self.mm2tom2


