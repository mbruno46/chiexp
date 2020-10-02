import numpy
import sys
import os

from scipy.optimize import minimize
 
p = os.path.realpath(__file__)
sys.path.append(f'{p[:-10]}/../../lib/python')
from chiexp import chisquare

[x0, x1, yy, dy] = numpy.loadtxt(f'{p[:-10]}/../averages.dat',unpack=True)
xx= numpy.c_[x0,x1]

W = numpy.diag(1./dy**2)

#f = lambda x, c, p: sum([(c==i+1)*p[i]*(1+p[3]*x[0]) for i in range(4)])
def f(x, c, a1, a2, a3, b):
    return locals()[f'a{int(c)}']*(1+b*x)

def df(x, c, a1, a2, a3, b):
    ai = []
    for i in [1,2,3]:
        ai.append(locals()[f'a{i}'])
    return [(c==i+1)*(1+b*x) for i in range(3)] + [sum([(c==i+1)*ai[i]*x for i in range(3)])]

cs = chisquare(xx,yy,W,f,df,v='x,c')

# runs the fit
cs.fit([1,1,1,1],minimize)

derfit = cs.derfit()

# we take to be uncorrelated
covmat = numpy.diag(dy**2)

# propagate errors of data y to parameters
covpars = derfit.T @ numpy.diag(dy**2) @ derfit
dpars= numpy.sqrt(numpy.diag(covpars))

try:
    import matplotlib.pyplot as plt

    plt.figure()
    xax = numpy.arange(0,max(x0),max(x0)/100)
    plt.title('Example of chiral & continuum fit')
    for i in range(3):
        idx=(x1==(i+1))
        plt.errorbar(x0[idx],yy[idx],dy[idx],fmt='.',color=f'C{i}')
        plt.plot(xax,[f(x,i+1,*cs.p) for x in xax],'-',color=f'C{i}')
    plt.show()
except:
    pass