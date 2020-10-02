# ChiExp

[![Python](https://img.shields.io/badge/Python-3.6-brightgreen.svg)](https://www.python.org)
[![Matlab](https://img.shields.io/badge/MATLAB-R2019b-brightgreen.svg)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

A package to compute the expectation value of the chi squared defined 
from arbitrary regularizations of the inverse covariance matrix, following the results of 
Ref. [[1][1]].

If you use this library in your publications please consider citing:

 - [[1][1]] M. Bruno, R. Sommer, In preparation.

### Authors

Copyright (C) 2017-20 Mattia Bruno, Rainer Sommer

## Installation

A MATLAB version of the package can be found in the 
directory [`/lib/matlab`](./lib/matlab); it can be imported
in any script by typing

```matlab
>> addpath '/path/to/chiexp/directory/lib/matlab'
>> help chiexp
```

We also provide a Python module contained
in the directory [`lib/python`](./lib/python), that can be easily 
imported with

```python
>>> import sys
>>> sys.path.append('/path/to/chiexp/directory/lib/python')
>>> from chiexp import chisquare
>>> help(chisquare)
```

## Documentation

Examples for can be found [here](./examples/), while
the documentation of the package can be found in [HTML][2]
or in [PDF][3] format

## Task list

 - [ ] support multiple ensembles python
 - [ ] support multiple ensembles matlab
 
[1]: https://arxiv.org
[2]: https://mbruno46.github.io/chiexp
[3]: ./docs/chiexp-doc.pdf
