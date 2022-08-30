import numpy
import sys
import os

from scipy.optimize import minimize
 
p = os.path.realpath(__file__)
sys.path.append(f'{p[:-10]}/../../lib/python')
from chiexp import chisquare

dat = numpy.loadtxt('../freefield-32x12-m0.25-r001.dat.gz')
(ncnfg, T) = dat.shape

xx=numpy.arange(T//2+1)
yy=numpy.mean(dat,axis=0)[0:T//2+1]
dy=numpy.sqrt(numpy.var(dat,axis=0))[0:T//2+1]

f = lambda x, a, m: a*(numpy.exp(-m*x) + numpy.exp(-m*(T-x)))
df = lambda x, a, m: [(numpy.exp(-m*x) + numpy.exp(-m*(T-x))),
                      a*(-m*numpy.exp(-m*x) + m*numpy.exp(-m*(T-x)))]

x0min = range(T//4)
p0 = [yy[0], 0.2]

res = []

for x0 in x0min:
    idx=slice(x0,T//2+1,None)
    print(f'\nFit with range [{x0}, {xx[-1]}]')
    
    W = numpy.diag(1./dy[idx]**2)
    
    cs = chisquare(xx[idx],yy[idx],W,f,df)
    cs.fit(p0, minimize)
   
    print(f'parameters = {cs.p}')
    print(f'chi^2 = {cs.c2:g}')
    
    [ce,dce,_] = cs.chiexp(dat[:,idx])
    print(f'chiexp         = {ce:g} +- {dce:g}')
    
    [p,dp,h] = cs.pvalue()
    print(f'chiexp from MC = {numpy.mean(h):g} +- {numpy.sqrt(numpy.var(h)/len(h)):g}')

    res.append([cs.c2, ce, dce, p, dp])

res = numpy.array(res)

try:
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.xlabel('x0min')
    plt.ylabel('Chi^2 / ChiExp')
    plt.errorbar(x0min,res[:,0]/res[:,1],res[:,0]/res[:,1]**2*res[:,2],fmt='.')
    plt.show()
    
    plt.figure()
    plt.xlabel('x0min')
    plt.ylabel('Quality of fit')
    plt.errorbar(x0min,res[:,3],res[:,4],fmt='.')
    plt.show()
except:
    pass