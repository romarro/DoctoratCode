# -*- coding: utf-8 -*-
"""
Created on Fri Jun 06 19:41:15 2014

@author: Vlad

luata de pe web:
http://blog.enthought.com/general/visualizing-uncertainty/#.U5Ht9PnCaSp
"""
from sklearn.gaussian_process import GaussianProcess
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

def pdense(x, y, sigma, M=1000):
    """ Plot probability density of y with known stddev sigma
    """
    assert len(x) == len(y) and len(x) == len(sigma)
    N = len(x)
    # TODO: better y ranging
    ymin, ymax = min(y - 2 * sigma), max(y + 2 * sigma)
    yy = np.linspace(ymin, ymax, M)
    a = [np.exp(-((Y - yy) / s) ** 2) / s for Y, s in zip(y, sigma)]
    A = np.array(a)
    A = A.reshape(N, M)
    plt.imshow(-A.T, cmap='gray', aspect='auto',
               origin='lower', extent=(min(x)[0], max(x)[0], ymin, ymax))
    plt.title('Density plot')

def gpr(seed=0, N=20, M=1000, sigma=1.0):
    """ from scikits.learn demo
    """
    np.random.seed(seed)

    def f(x):
        """The function to predict."""
        return x * np.sin(x)
    
    X = np.linspace(0.1, 9.9, 20)
    X = np.atleast_2d(X).T
    y = f(X).ravel()
    y = np.random.normal(y, sigma)
    x = np.atleast_2d(np.linspace(0, 10, M)).T
    nugget = (sigma / y) ** 2
    gp = GaussianProcess(corr='squared_exponential', theta0=1e-1,
                         thetaL=1e-1, thetaU=1.0,
                         nugget=nugget,
                         random_start=100)
    gp.fit(X, y)
    y2, MSE = gp.predict(x, eval_MSE=True)
    s2 = np.sqrt(MSE)
    return X, y, x, y2, s2
if __name__=='__main__':
    
    X, y, x, y2, s2 = gpr(seed=0)
    plt.figure(1)
    pdense(x, y2, s2, M=1000)
    plt.plot(X, y, 'r.')
    plt.plot(x, y2, 'b:')
    a = plt.gca()
    a.set_ylim(-10, 15)
    plt.xlabel('$x$')
    plt.ylabel('$f(x)$')
    plt.show()