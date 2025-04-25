.. highlight:: python
   :linenothreshold: 3

Python documentation
====================

The package `chiexp` can be loaded in python using the
standard import command::

   import sys
   sys.path.append('/path/to/chiexp/directory/lib/python')
   from chiexp import chisquare


The chisquare class
-------------------

.. automodule:: chiexp
   :members: chisquare

Examples of usage
-----------------

First we define a new instance of the `chisquare` class; in the 
following example we do it for an uncorrelated fit

.. code-block:: python

   [x, y, dy] = load_your_data_set() 
   
   # if dy is the error of y we define W for an uncorrelated fit
   W = [1./e**2 for e in dy]
   
   func = lambda x, a, m: a*exp(-m*x)
   dfunc = lambda x, a, m: [exp(-m*x), -a*x*exp(-m*x)]

   c=chisquare(x,y,W,func,dfunc)


If the minimum of :math:`\chi^2` is not know the user can use the 
``fit`` method; otherwise the values of 
the parameters at the minimum can be passed to the class
via the ``chisq`` method

.. code-block:: python

  >>> from scipy.optimize import minimize
  >>> p0=[1.,1.] # guess values of the paramters
  >>> [p, c2] = c.fit(p0, minimize)
  >>> # or pass the parameters at the minimum
  >>> c.chisq(p)


The user can also compute the error of the fitted 
parameters from the errors (fluctuations) of the input 
observables :math:`Y`. To do so the derivatives 
:math:`d p_\alpha / d y_i` (defined at the minimum of 
the :math:`\chi^2`) are needed and if the covariance
matrix is known, applying the chain rule returns the wanted
errors

.. code-block:: python

   >>> der = c.derfit() # N x Na matrix
   >>> print (der.T @ cov @ der) # errors of parameters


:math:`\langle \chi^2 \rangle` can be computed using the 
method `chiexp`

.. code-block:: python

  >>> [ce, dce, _] = c.chiexp(cov)
  >>> # if cov contains the fluctuations of M configs and 2 replica
  >>> (M, N) = numpy.shape(cov)
  >>> nr = [N1 N2] # such that M=N1+N2
  >>> print M - sum(nr)
  0
  >>> [ce,dce,covest] = c.chiexp(cov,Nrep=nr,Stau=2.5)


.. The input argument `cov` can be either a matrix (`list`
   or `numpy.ndarray`) of dimensions NxN or a matrix 
   of dimensions MxN. In the first case the program 
   assumes that `cov` corresponds to the covariance matrix
   previously estimated by the user. In the second case
   the program assumes that `cov` contains the fluctuations
   of the observable `y` (mean value subtracted from each
   measurement) over M different configurations: here 
   the programs computes a *derived autocorrelation function*
   that is used to estimate :math:`\langle\chi^2\rangle` and its error.

Finally the quality of fit is estimated from a Monte-Carlo
chain of length `nmc`

.. code-block:: python

  >>> nmc=10000
  >>> [p, dp, h] = c.pvalue(nmc, plot=True)


The method `pvalue` returns the quality of fit and its error,
`p` and `dp` respectively, and the Monte Carlo history of 
:math:`\chi^2` `h`. If the `plot` flag is activated a plot is automatically
generated with the distribution probability, namely a
normalized histogram of `h`.

