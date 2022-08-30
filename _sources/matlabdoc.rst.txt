.. highlight:: matlab
   :linenothreshold: 3

Matlab documentation
====================

The matlab package contains three functions:

:derfit: computes the derivatives of the parameters w.r.t. the input observable Y at the minimum

:chiexp: computes the expected :math:`\chi^2`

:qfit: computes the quality of fit


derfit
------

At the minimum of the :math:`\chi^2`, the parameters depend on the input 
observables Y. This function computes their derivatives w.r.t. Y::

   der = derfit(X,Y,W,p,'field',value);


:Parameters: 
  - **X, Y** - arrays with N entries; Y is assumed to contain the 
    central values of the observable
  
  - **W** - weight matrix, N-by-N, used in the fit
  
  - **p** - array with the values of the parameters at the minimum

  - **field** - a string whose admitted values as ``'f'`` and ``'df'``

  - **value** - the value corresponding to one of the two strings above

     - **df** - the gradient of the function f. An array of functions is expected 
       ``df = @(x,p) [df/dp1 df/dp2 ...]``
     
     - **f** - fitted function. It must be passed in the form ``f=@(x,p) f(x,p)``;
       if passed together with `df` it is used to check numerically the 
       provided gradient; if passed without `df` it is used to numerically compute 
       the gradient

:Additional arguments:
  
  - **c** *(optional)*: array of length N used as additional argument for 
    `f` and `df`, which 
    in this case must obey the syntax ``f = @(x,p,c) f(x,p,c)`` (same for `df`);
    it can be useful for global fits

:Returns:
  - **der**: a N-by-NA matrix, with NA the number of parameters; contains the derivatives
    of the parameters ``der(i,a) = dp(a)/dy(i)``


**Examples**::

   % derfit can be used in two ways:
   % either by passing the function (numerical gradient)
   func = @(x,p) p(1) + p(2)*x;
   der = derfit(X,Y,W,p,'f',func);
   
   % or by passing the gradient (analytical)
   grad = @(x,p) [1, x];
   der = derfit(X,Y,W,p,'df',grad);
   
   % if both are passed, the function is used to check the
   % gradient
   der = derfit(X,Y,W,p,'f',func,'df',grad);

   % additional argumemts
   der = derfit(X,Y,W,p,'f',func,'df',grad,'c',c);


chiexp
------

This function computes the expected :math:`\chi^2`. The covariance matrix is 
automatically estimated from the fluctuations of Y, but the user can also
provide it independently::

   [ce,dce,nu,covest] = chiexp(X,Y,W,p,cov,'field',value)


:Parameters:
  - **X, Y** - arrays with N entries
  
  - **W** - weight matrix. It can be an array of dimensions N or a matrix of 
    dimension NxN. In the first case a diagonal weight matrix is automatically
    created
  
  - **p** - array with the values of the parameters at the minimum

  - **cov** - if cov is a N-by-N matrix, the programs takes it as the covariance 
    matrix, assuming the user has previously computed it; if cov is a
    M-by-N matrix, the program assumes that it contains the values of the
    N observables Y on M different configurations.

  - **field** - a string whose admitted values as ``'f'`` and ``'df'``

  - **value** - the value corresponding to one of the two strings above

      - **f** - fitted function. It must be passed in the form ``f = @(x,p) f(x,p)``; 
        this argument is ignored if `df` is passed
      
      - **df**: the gradient of the function f. An array of functions is expected
        ``df = @(x,p) [df/dp1 df/dp2 ...]``. If empty, the gradient is computed 
        numerically and *f* must be passed


Additional parameters can be passed using the syntax::
   
   [ce,dce] = chiexp(X,Y,W,p,cov,'df',df,'field1',value1,'field2',value2...)
   % or 
   [ce,dce] = chiexp(X,Y,W,p,cov,'f',f,'field1',value1,'field2',value2...)


:Additional arguments: 
  - **c** *(optional)*: array of length N used as additional argument for 
    `f` and `df`, which 
    in this case must obey the syntax ``f = @(x,p,c) f(x,p,c)`` (same for `df`);
    it can be useful for global fits
  
  - **Nrep** *(optional)*: array with the number of configurations belonging to
    each replica, e.g. ``[N1 N2 N3 ...]``. 
    If not given, the number of replica is assumed to be 1.
  
  - **Stau** *(optional)*: value of the parameter Stau. Default is 1.5.
  
  - **plot** *(optional)*: flag to produce a plot of the expected 
    :math:`\chi^2` as a function of the
    window. Default is ``'off'``. Accepted values are ``'on'`` and ``'off'``.

:Returns:
  - **ce, dce** - the values of :math:`\langle \chi^2 \rangle` and its error

  - **nu**: the matrix :math:`\big[ C^{1/2} W^{1/2} (1-P) W^{1/2} C^{1/2} \big]`

  - **covest**: the estimated covariance matrix using the window obtained 
    from the expected :math:`\chi^2` (which may not be adeguate for general purposes)

**Examples**::
   
   % chiexp can be used in two ways:
   % either by passing the function (numerical gradient)
   func = @(x,p) p(1) + p(2)*x;
   [ce,dce] = chiexp(X,Y,W,p,cov,'f',func);
   
   % or by passing the gradient (analytical)
   grad = @(x,p) [1, x];
   [ce,dce] = chiexp(X,Y,W,p,cov,'df',grad);

   % if cov is M-by-N and has 2 replica
   Nrep=[N1, N2]; % M=N1+N2
   [ce,dce] = chiexp(X,Y,W,p,cov,'df',grad,'Nrep',Nrep,'Stau',2.5);

   % full output
   [ce,dce,nu,covest] = chiexp(X,Y,W,p,cov,'df',grad);

   % for W^-1 = diag(C)
   [ce,dce,nu,covest,ce2,dce2] = chiexp(X,Y,W,f,p)


.. note:: 
   If `W` is an array of length N corresponding to the diagonal
   entries of the covariance matrix, the expected :math:`\chi^2` can be obtained with a 
   second formula. To have both results the user must request two additional 
   output arguments, that will contain 
   :math:`\langle \chi^2 \rangle = N - \mathrm{tr}\big[C W^{1/2} P W^{1/2}\big]`
   and its error 


qfit
----

The quality of fit and its error are estimated from a Monte Carlo integration.::

   [Q,dQ] = qfit(nmc,nu,c2,plot)

**Input arguments**:

* `nmc`: the size of the Monte Carlo chain

* `nu`: the matrix returned from the function `chiexp`

* `c2`: the value of the :math:`\chi^2` obtained from the fitting procedure

* `plot` (*optional*): is passed a histogram with normalized distribution
   probability obtained from the Monte Carlo process is plotted

**Returns**:

* `Q,dQ`: the quality of fit and its error
