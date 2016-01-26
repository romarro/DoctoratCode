# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 14:24:40 2014

@author: Vlad
"""

import numpy as np
import uncertainties as un
import uncertainties.numpy as unp
import csv


x=np.array(un.ufloat(1,0.01),un.ufloat(2,0.02),un.ufloat(3.1,0.01))

writer=csv.writer()