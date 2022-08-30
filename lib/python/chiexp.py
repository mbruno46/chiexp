#################################################################################
#
# chiexp.py
# Copyright (C) 2017-22 Mattia Bruno, Rainer Sommer
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#################################################################################

import numpy
import inspect

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB=True
except:
    MATPLOTLIB=False

__all__=['chisquare']


class ChiExpError(Exception):
    pass

def _gamma(dat,Nrep):
    # N is the number of configs
    (N,M)=numpy.shape(dat)
    Wmax=min(Nrep)//2
    g = numpy.zeros((M,M,Wmax+1), dtype=numpy.double)
    ofs=0
    for Nr in Nrep:
        imax=ofs+Nr
        aux=numpy.fft.fft(dat[ofs:imax,:],n=2*N,axis=0)
        aux2=numpy.zeros((2*N,M,M),dtype=numpy.float64) # note that this forces line below to pick only real part
        for i in range(2*N):
            aux2[i,:,:] = numpy.outer(aux[i,:],numpy.conj(aux[i,:])).real
        aux2 = numpy.fft.ifft(aux2, axis=0)
        for a in range(M):
            for b in range(M):
                g[a,b,:] += aux2[0:Wmax+1,a,b].real
        ofs+=Nr
    n = 1./(N-len(Nrep)*numpy.arange(0.,Wmax+1,1.))
    for i in range(M):
        for j in range(i,M):
            g[i,j,:] *= n
            if i!=j:
                g[j,i,:] = g[i,j,:]
    return g


def eval_chiexp(gam,mat,ncnfg):
    N=numpy.shape(gam)[2]
    c=numpy.zeros((N,))
    for i in range(N):
        c[i] = numpy.trace(mat @ gam[:,:,i])
    return c/(float)(ncnfg)
    
    
def eval_chiexp_fast(dat, Nrep, mat):
    w,v = numpy.linalg.eig(mat)

    # N is the number of configs
    (N,M)=numpy.shape(dat)
    Wmax=min(Nrep)//2
    g = numpy.zeros((M,Wmax+1), dtype=numpy.double)
    ofs=0
    for Nr in Nrep:
        imax=ofs+Nr
        dy = dat[ofs:imax,:]
        
        aux = numpy.fft.fft(dy @ v, n=2*Nr, axis=0)
        aux2 = numpy.fft.ifft(aux * numpy.conj(aux), axis=0)
        for a in range(M):
            g[a,:] += aux2[0:Wmax+1,a].real
        ofs += Nr
        
    n = 1./(N-len(Nrep)*numpy.arange(0.,Wmax+1,1.))
    return (w @ g).real * n / N

    
def _find_window(rho, N, Stau):
    Wmax = int(len(rho))
    rho_int = 0.
    flag=0
    for W in range(1,Wmax):
        rho_int += rho[W]
        tauW = Stau/numpy.log(numpy.fabs((rho_int+1.)/rho_int))
        gW = numpy.exp(-W/tauW) - tauW/numpy.sqrt(W*N)
        if (gW<0):
            Wopt = W
            Wmax = min([Wmax,2*Wopt])
            flag=1
            break
    if (flag==0):
        print('Warning: automatic window procedure failed')
        Wopt = Wmax
    else:
        print(f'Automatic window set at {Wopt}')
    return [Wopt, Wmax]


class chisquare:
    """ A class with the relevant functionalities to compute the 
    expected chi square.

    Parameters:
       x (array): array with the values of the x-axis; 2-D arrays are
          accepted and the second dimension is intepreted as internal kinematic
          index
       y (array): array with the values of the y-axis, of same length as `x`
       W (array): the weight matrix, N-by-N, or its diagional of length N
       f (function): callable function or lambda function defining the fitted function;
          the program assumes `x` correspond to the first arguments
       df (function): callable function or lambda function returning an array
          that contains the gradient of `f`, namely :math:`\partial
          \phi(\{p\},\{x_i\})/\partial p_\\alpha`
       v (str, optional): a string with the list of variables used in `f` as
          the kinematic coordinates. Default value corresponds to `x`, which
          implies that `f` must be defined using `x` as first and unique
          kinematic variable.
          
    """

    def __init__(self,x,y,W,f,df,v='x'):
        if not isinstance(x,numpy.ndarray):
            raise ValueError('x must be a numpy.array')
        if not isinstance(y,numpy.ndarray):
            raise ValueError('y must be a numpy.array')

        if numpy.ndim(y)>1:
            raise ValueError('y should be a 1D array')
            
        if numpy.ndim(x)==1:
            self.n=len(x)
            self.nx=1
        elif numpy.ndim(x)==2:
            (self.n, self.nx) = numpy.shape(x)
        else:
            raise ValueError('x should be a 1D or 2D array')
            
        self.x = numpy.reshape(x,(self.n,self.nx))
        if len(y)!=self.n:
            raise ValueError('x and y lengths do not match')
        self.y = numpy.array(y)
        
        if isinstance(W,(list,numpy.ndarray)):
            if numpy.ndim(W)==1:
                self.W = numpy.diag(W)
            else:
                self.W = numpy.array(W)

        # checks W is positive definite
        eig=numpy.linalg.eig(self.W)[0]
        if numpy.any(eig<0)==True:
            print('WARNING: W is not a positive definite matrix')
        
        self.v = v.rsplit(',')
        self.f = f
        self.df = df
        args = inspect.getargspec(f)[0]
        if args != inspect.getargspec(df)[0]:
            raise ChiExpError(f'Unexpected f and df: arguments do not match')
        self.pars = []
        for vn in args:
            if not vn in self.v:
                self.pars.append(vn)
        self.np = len(self.pars)
        self.e = numpy.zeros((self.n,),dtype=numpy.float64)
        self.p = [0.0]*self.np
        self.c2=None


    def chisq(self,p):
        """ 
        Returns the chi square from the input parameters `p`.
        """
        self.set_pars(p)
        for i in range(self.n):
            self.e[i] = self.f(*self.x[i,:], *self.p) - self.y[i]
        return self.e @ self.W @ self.e
   

    def fit(self,p0,min_search):
        """ 
        Minimizes the chi square starting from the initial guess parameters `p0`
        
        Parameters:
            p0 (array): the initial guess values of the parameters
            min_search (function): an external minimizer
            
        Returns:
            list: output parameters and the value of the chi square at the minimum
        """
        res = min_search(lambda p:self.chisq(p),p0)
        self.c2 = res.fun
        self.set_pars(res.x)
        return [res.x, res.fun]


    def set_pars(self,p0):
        """
        Sets the parameters to the input values `p0`
        """
        if len(p0)!=self.np:
            raise ChiExpError('unexpected number of parameters')
        for i in range(self.np):
            self.p[i] = p0[i]


    def derfit(self):
        """ 
        Returns the derivative of the parameters w.r.t. Y

        At the minimum the parameters depend on the input observables; their
        errors can be obtained by error propagation through the derivatives
        provided by this routine.

        Returns:
            array: the derivatives ``der(i,a) = dp(a)/dy(i)``
        """
        g=numpy.array([self.df(*self.x[i,:], *self.p) for i in range(self.n)]) # N x Na matrix
        Hmat = g.T @ self.W @ g # Na x Na matrix
        Hinv=numpy.linalg.inv(Hmat) 
        return self.W @ g @ Hinv

    
    def chiexp(self,cov,Nrep=None,Stau=1.5,Wopt=None,Wcov=None,plot=False):
        """
        Computes the expected chi square

        Parameters:
           cov (list or array): the covariance matrix is assumed if the object passed has 
                                is a N-by-N 2D array; otherwise if it is a M-by-N 2D array the program
                                assumes that it contains the fluctuations of the N observables
                                over M configurations
           Nrep (list, optional): number of configurations per replica, e.g. if ``Nrep=[N1,N2,N3]`` then
                                  ``M=N1+N2+N3``
           Stau (float, optional): parameter used in the automatic window procedure
           Wopt (float, optional): optimal window used in the estimate of expected chi square; if passed `Stau` is ignored
           Wcov (float, optional): window used in the estimate of all entries of the covariance matrix
           plot (bool, optional): if set to `True` the program produces a plot with the autocorrelation 
                                  function of the expected chi square

        Returns:
           list: the expected chi square, its error and the estimated covariance matrix; if `Wcov` is not 
                 passed the window obtained from the autocorrelation function of the expected chi square is used 
                 for all elements; if `cov` is a N-by-N 2D array it returns a copy of it

        Note:
           The additional parameters `Nrep`,`Stau` and `plot` are ignored if `cov` 
           corresponds to the covariance matrix, i.e. it is a N-by-N array.
           The matrix :math:`[C^{1/2} W^{1/2} (1-P) W^{1/2} C^{1/2}]` can be accessed
           by invoking the ``nu`` element of the ``chisquare`` class
        """
        g=numpy.array([self.df(*self.x[i,:], *self.p) for i in range(self.n)]) # N x Na matrix
        Wg=self.W @ g

        Hmat = g.T @ self.W @ g # Na x Na matrix
        Hinv=numpy.linalg.inv(Hmat) 

        PP=self.W - Wg @ Hinv @ Wg.T
    
        if isinstance(cov,(list,numpy.ndarray)):
            if numpy.shape(cov)==(self.n,self.n):
                _cov=numpy.array(cov)
                ce=numpy.trace(PP.dot(_cov))
                dce=0
                wopt=None
            elif numpy.shape(cov)[1]==self.n:
                ncnfg=numpy.shape(cov)[0]
                _yy=numpy.mean(cov,axis=0)
                _y = numpy.array(cov) - numpy.array([_yy]*ncnfg) # computes fluctuations
                
                if Nrep is None:
                    Nrep=[ncnfg]
                gg=_gamma(_y,Nrep)
#                 chiexp_t = eval_chiexp(gg,PP,ncnfg)
#                 print(numpy.sum(chiexp_t - eval_chiexp_fast(_y,Nrep,PP)))
                chiexp_t = eval_chiexp_fast(_y,Nrep,PP) 
                
                if not Wopt is None:
                    wopt=Wopt
                    wmax=2*Wopt
                else:
                    # uses normalized autocorr function
                    [wopt, wmax]=_find_window(chiexp_t/chiexp_t[0],ncnfg,Stau)
                ce=chiexp_t[0]+2.*numpy.sum(chiexp_t[1:wopt+1])
                dce=numpy.sqrt((4*wopt+2)/ncnfg)*ce
                
                if Wcov is None:
                    Wcov=wopt

                _cov=gg[:,:,0]+2.*numpy.sum(gg[:,:,1:Wcov+1],axis=2)
                _cov /= (float)(ncnfg)

                if plot and MATPLOTLIB:
                    rho=numpy.r_[chiexp_t[0:wmax]/chiexp_t[0],[0.]*(wmax+wopt)]
                    drho=numpy.zeros((wmax,))
                    for t in range(wmax):
                        for k in range(max(1,t-wopt),t+wopt):
                            drho[t]+=(rho[k+t+1]+rho[abs(k-t)+1]-2*rho[t+1]*rho[k+1])**2;
                    drho=numpy.sqrt(drho/(float)(ncnfg))

                    plt.figure()
                    plt.title('ChiExp automatic window procedure')
                    plt.xlabel('[config index]')
                    plt.plot([0,wmax],[0,0],'-k',lw=.75)
                    plt.errorbar(range(wmax),chiexp_t[0:wmax]/chiexp_t[0],drho,fmt='-o',label='Gamma(t)/Gamma(0)')
                    plt.plot([wopt,wopt],[0,0.8],'-r',label='Wopt')
                    plt.xlim(0,wmax)
                    plt.legend(loc='upper right')
                    plt.show()
            else:
                raise ValueError('Unexpected cov, neither %dx%d nor Mx%d matrix' % (self.n,self.n,self.n))
        else:
            raise ChiExpError('Unexpected cov, neither %dx%d nor Mx%d matrix' % (self.n,self.n,self.n))

        # checks cov matrix before taking square root
        [w,v] = numpy.linalg.eig(_cov)
        if numpy.any(w<0):
            print('The estimated covariance matrix has negative eigenvalues with automatic window = %d' % (Wcov))
            mask = w>1e-12
            Csqroot = v[:,mask] @ numpy.diag(numpy.sqrt(w[mask])) @ v[:,mask].T
        else:
            Csqroot= v @ numpy.diag(numpy.sqrt(w)) @ v.T

        self.nu= Csqroot @ PP @ Csqroot
        self.ce = [ce,dce]

        self.PP = PP
        self.da = Hinv @ Wg.T
        return [ce,dce,_cov]

    def dchiexp(self):
        return numpy.sqrt(2.*numpy.trace(self.nu @ self.nu))
    
    def pvalue(self,method='eig',nmc=5000,plot=False):
        """ 
        Computes the p-value of the fit

        Parameters:
           method (string, optional) : string specifying the method to estimate the 
              quality of fit. Accepted values are 'MC' for a pure Monte Carlo estimate 
              or 'eig' (default) for the formula based on the eigenvalues of the matrix `nu`. 
           nmc (int, optional) : number of Monte Carlo samples used to estimate the quality of fit.
              Default is 5000.
           plot (bool, optional): if set to True plots the probabilty distribution of the 
              expected chi square

        Returns:
           list: the quality of fit, the error of the quality of fit and the Monte Carlo history of the expected chi square

        Note:
           The class must know the value of the chi square at the
           minimum, which means that either `fit` or `chisq` must 
           be called before `qfit`.
           If method is set to 'MC' the error of the quality of fit is based only the 
           MC sampling. 
        """
        cexp=[]
        dcexp=[]

        if method=='MC':
            z=numpy.random.normal(0.,1.,self.n*nmc).reshape(nmc,self.n)
            cexp = numpy.einsum('ia,ab,bi->i',z,self.nu,z.T)
        elif method=='eig':
            ev=numpy.linalg.eig(self.nu)[0]
            ev=ev.real
            
            # sorts eigenvalues, such that first is positive and largest
            ev=numpy.sort(ev)[::-1]
            
            # chops eigenvalues whose real part is smaller than eps
            eps=1e-14 * numpy.max(numpy.abs(ev))
            ev=ev*(ev>eps)
            
            for i in range(nmc):
                z=numpy.random.normal(0.,1.,self.n)
            
                cexp.append( ev @ z**2 )
                
                x = (self.c2 - ev[1:] @ z[1:]**2)/ev[0]
                if x>0:
                    dcexp.append( numpy.exp(-x*0.5)*0.5/numpy.sqrt(x) )
                else:
                    dcexp.append(0.)
            
        th=numpy.array(cexp)<self.c2
        p=1.0 - numpy.mean(th)
        dp=numpy.std(th,ddof=1)/(nmc)**0.5
        
        if plot and MATPLOTLIB:
            plt.figure()
            plt.title('Probability distribution')
            plt.ylabel('$P(\chi^2)$')
            plt.xlabel('$\chi^2$')
            h = plt.hist(cexp,density=True,bins=40,label='MC')
            plt.plot([self.c2,self.c2],[0,max(h[0])],label='$\chi^2$')
            plt.legend(loc='upper right')
            plt.show()

        return [p,dp,numpy.array(cexp)]
