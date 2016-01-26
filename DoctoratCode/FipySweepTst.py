from __future__ import division
from fipy import *
import numpy as np

"""
Rezolvam problema difuziei in cazul in care coeficientul de difuzie depinde de variabila.

\frac{\partial \phi}{\partial t}=\nabla [ D(\phi) \cdot \nabla \phi ]

D=D0*(1-\phi)

"""

nx = 50
dx=1.
mesh=Grid1D(dx=dx,nx=nx)
valueLeft = 1.
valueRight =0.
steps=100
D0=1.
timestepduration =0.9*dx**2/(2*D0)

phi = CellVariable(name='$\phi$',mesh=mesh,value=valueRight,hasOld=1)
eq = TransientTerm() == DiffusionTerm(D0*(1-phi))
phi.constrain(valueLeft,mesh.facesLeft)
phi.constrain(valueRight,mesh.facesRight)

if __name__ == '__main__':
    viewer = Viewer(phi,datamin=0.,datamax=1.)
    viewer.plot()

for sweep in range(5):
    raw_input("step")
    rez = eq.sweep(var=phi,dt=timestepduration)
    if __name__=='__main__':
        viewer.plot()
    print rez

raw_input("gata")
