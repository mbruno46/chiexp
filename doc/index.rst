ChiExp Documentation
====================

A package to compute the expectation value of the :math:`\chi^2` defined 
from arbitrary weight matrices :math:`W`

.. math::
   \chi^2 = \sum_{i,j} [y_i - f(x_i,a)] W_{ij} [y_j - f(x_j,a)]


The implementation is based on the results of Ref. [1]_
and we summarize below the main equation

.. math::
   \langle \chi^2 \rangle = \mathrm{tr} \ \big[C_W W^{1/2} (1-P) W^{1/2} \big] \,,
   \quad C_W = W^{1/2} C W^{1/2}

with :math:`C` being the covariance matrix and :math:`P` a projector, 
defined from the function :math:`f`
and its derivatives computed at the minimum of the :math:`\chi^2`
(for more details read Ref. [1]_).

Autocorrelations can be taken into account in a straight-forward manner by replacing

.. math::
   C_{ij} = \frac{1}{N} \sum_{t=-\infty}^\infty \Gamma_{ij}(t) \,,

with :math:`N` the number of configurations and :math:`\Gamma` the autocorrelation function.


The package is written primarily for Matlab, but we provide also a Python
implementation in the directory `chiexp`.
Examples for Matlab can be found in the folder `examples`.

..
.. toctree::
   :maxdepth: 1
   :caption: The documentation can be found here

   matlabdoc
   pydoc


Further details can be found by typing `help chiexp` in Matlab
or `help(chisquare)` in Python.

References
----------

.. [1] M. Bruno and R. Sommer `title <https://arxiv.org>`__

Differences
-----------

1. derivatives :math:`df(x_i,a)/da_\alpha`: in the Python implementation the derivatives
   are computed analytically using the `sympy` package; in the Matlab library
   the user can either provide them (analytically), or obtain the numerical estimates
   by passing the function instead. If they are passed directly the library can also
   perform a numerical check.

2. in Python a unique class `chisquare` is initialized by the user, which contains all
   the important methods. In Matlab several functions are provided and the user must
   follow a simple workflow.
