import numpy
import sys
import os

p = os.path.realpath(__file__)
sys.path.append(f'{p[:-10]}/../../lib/python')
from chiexp import chisquare

am=0.5
T=16
ncnfg=1000

# analytical prediction with infinite Time extent
omega = numpy.arccosh(0.5*am**2+1.)
Gt = lambda x: numpy.exp(-abs(x)*omega)/(2*numpy.sinh(omega))

# covariance matrix
cov_theo = numpy.zeros((T,T));
for t1 in range(T):
    for t2 in range(T):
        cov_theo[t1,t2] = Gt(t1)*Gt(t2) + Gt(t1-t2)*Gt(0)

print('eigenvalues of exact cov matrix')
print(numpy.linalg.eig(cov_theo)[0])

xx = numpy.arange(T)
yy = numpy.array([Gt(x) for x in xx])
dy = numpy.sqrt(numpy.diag(cov_theo))

f = lambda x, a, m: a*numpy.exp(-m*x)
df = lambda x, a, m: [numpy.exp(-m*x), -x*a*numpy.exp(-m*x)]

p0 = [1./(2.*numpy.sinh(omega)), omega]
print(f'amp  = {p0[0]}')
print(f'mass = {p0[1]}');

#####

print('\nW = inverse of diag of cov matrix')
W = numpy.diag(cov_theo)

cs = chisquare(xx,yy,W,f,df)
print(f'chi2 = {cs.chisq(p0)} / dof = {T-len(p0)}')

[ce, dce,_] = cs.chiexp(cov_theo)
print(f'chiexp = {ce:.4f} +- {dce:.4f}')
assert ce < T-len(p0)

#####

print('\nW = inverse of cov matrix')
W = numpy.linalg.inv(cov_theo)

cs = chisquare(xx,yy,W,f,df)
print(f'chi2 = {cs.chisq(p0)} / dof = {T-len(p0)}')

[ce, dce,_] = cs.chiexp(cov_theo)
print(f'chiexp = {ce:.4f} +- {dce:.4f}')
assert abs(ce - (T-len(p0))) < 1e-7