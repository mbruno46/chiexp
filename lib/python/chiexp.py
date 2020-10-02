####################################################
#
# Copyright (c) 2017-19 Mattia Bruno, Rainer Sommer.
#
####################################################

import numpy
import inspect
#from scipy.optimize import minimize
#from scipy.linalg import sqrtm
#import sympy
#from sympy.parsing.sympy_parser import parse_expr
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
        c[i] = numpy.trace(mat.dot(gam[:,:,i]))
    return c/(float)(ncnfg)


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
          
    Examples:
    >>> 
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
            list: output parameters
            float: value of the chi square at the minimum
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

    
    def chiexp(self,cov,Nrep=None,Stau=1.5,Wopt=None,plot=False):
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
            plot (bool, optional): if set to `True` the program produces a plot with the autocorrelation 
                                   function of the expected chi square

        Returns:
            float: the expected chi square
            float: the error of the expected chi square
            array: if cov contains the fluctuations, this is the estimated 
                   covariance matrix using the same window for all elements 
                   from the expected chi square; otherwise it is a copy of `cov`

        Note:
            The additional parameters `Nrep`,`Stau` and `plot` are ignored if `cov` 
            corresponds to the covariance matrix, i.e. it is a N-by-N array.
        
        Note:
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
            elif numpy.shape(cov)[1]==self.n:
                ncnfg=numpy.shape(cov)[0]
                _yy=numpy.mean(cov,axis=0)
                _y = numpy.array(cov) - numpy.array([_yy]*ncnfg) # computes fluctuations
                
                if Nrep is None:
                    Nrep=[ncnfg]
                gg=_gamma(_y,Nrep)
                chiexp_t = eval_chiexp(gg,PP,ncnfg)
                if not Wopt is None:
                    wopt=Wopt
                    wmax=2*Wopt
                else:
                    # uses normalized autocorr function
                    [wopt, wmax]=_find_window(chiexp_t/chiexp_t[0],ncnfg,Stau)
                ce=chiexp_t[0]+2.*numpy.sum(chiexp_t[1:wopt+1])
                dce=numpy.sqrt((4*wopt+2)/ncnfg)*ce

                _cov=gg[:,:,0]+2.*numpy.sum(gg[:,:,1:wopt+1],axis=2)
                _cov /= (float)(ncnfg)

                if plot:
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
                    ##plt.plot(range(len(chiexp_t)),chiexp_t/chiexp_t[0],'-ob',label='Gamma(t)/Gamma(0)')
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
            msg=f"""The estimated covariance matrix has negative eigenvalues with (automatic) window = {wopt}
                    Consider reducing the window by passing the parameter W"""
            raise ChiExpError('The estimate covariance matrix has negative eigenvalues with automatic window = %d' % (wopt))
        else:
            Csqroot= v @ numpy.diag(numpy.sqrt(w)) @ v.T
                    
        self.nu= Csqroot @ PP @ Csqroot
        self.ce = [ce,dce]
        
        return [ce,dce,_cov]

    def dchiexp(self):
        return numpy.sqrt(2.*numpy.trace(self.nu @ self.nu))
    
    def qfit(self,method='MC',nmc=5000,plot=False):
        """ Computes the quality of fit

        Parameters
        ----------
        method: string, optional
           string specifying the method to estimate the quality of fit. 
           Accepted values are 'MC' (default) for a pure Monte Carlo estimate 
           or 'eig' for the formula based on the eigenvalues of the matrix nu. 
        nmc: int, optional
           number of Monte Carlo samples used to estimate the quality of fit
        plot: bool, optional
           flag to plot the probabilty distribution of the expected chi square

        Returns
        -------
        float
            the quality of fit
        float
            its error
        numpy.ndarray
            the measurements of the expected chi square along the Monte Carlo history

        Note
        ----
        The class must know the value of the chi square at the
        minimum, which means that either `chiexp.fit` or `chiexp.pars` must 
        be called before qfit.

        """

        if self.pars is None:
            raise StandardError('the methods fit or pars should be called before qfit')

        cexp=[]
        cexp_i=[]
        dcexp=[]
        
        # checks if nu is potentially complex
        if numpy.trace(self.nu).imag>1e-13:
            raise ValueError('nu is complex')
        if numpy.trace(self.nu).real<0:
            raise ValueError('nu is negative')
            
        if method=='MC':
            for i in range(nmc):
                z=numpy.random.normal(0.,1.,self.n)
                cexp.append( z.dot(self.nu).dot(z).real )
                cexp_i.append( z.dot(self.nu).dot(z).imag )
        elif method=='eig':
            ev=numpy.linalg.eig(self.nu)[0]
            ev=ev.real
            evi=ev.imag
            
            # sorts eigenvalues, such that first is positive and largest
            ev=numpy.sort(ev)[::-1]
            
            # chops eigenvalues whose real part is smaller than eps
            eps=1e-14 * numpy.max(numpy.abs(ev))
            ev=ev*(ev>eps)
            
            for i in range(nmc):
                z=numpy.random.normal(0.,1.,self.n)
            
                cexp.append( ev.dot(z**2) )
                cexp_i.append( evi.dot(z**2) )
                
                x = (self.c2 - ev[1:].dot(z[1:]**2))/ev[0]
                if x>0:
                    dcexp.append( numpy.exp(-x*0.5)*0.5/numpy.sqrt(x) )
                else:
                    dcexp.append(0.)

        # mean of imag part < 2 * its error
        vi = numpy.abs(numpy.mean(cexp_i))
        ei = numpy.sqrt(numpy.var(cexp_i)/nmc)
        if vi>2*ei:
            raise ValueError('tr(nu) is complex, imaginary part = %.2e +- %.2e' % (vi,ei))
            
        th=numpy.array(cexp)<self.c2
        Q=1.0 - numpy.mean(th)
        dQ0=numpy.var(th)/nmc
        if method=='eig':
            dQ1=self.ce[1]/self.ce[0] * numpy.mean(dcexp)
        else:
            dQ1=0.0
        dQ=numpy.sqrt(dQ0 + dQ1**2)
        #if method=='MC':
        #    dQ=numpy.sqrt(numpy.var(th)/nmc)
        #elif method=='eig':
        #    dQ=self.ce[1]/self.ce[0] * numpy.mean(dcexp)
        
        if plot:
            plt.figure()
            plt.title('Probability distribution')
            plt.ylabel('$P(\chi^2)$')
            plt.xlabel('$\chi^2$')
            h = plt.hist(cexp,density=True,bins=40,label='MC')
            plt.plot([self.c2,self.c2],[0,max(h[0])],label='$\chi^2$')
            plt.legend(loc='upper right')
            plt.show()

        return [Q,dQ,numpy.array(cexp)]
