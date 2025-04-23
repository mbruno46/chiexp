# ChiExp

[![Test Octave](https://github.com/mbruno46/chiexp/actions/workflows/octave.yml/badge.svg)](https://github.com/mbruno46/chiexp/actions/workflows/octave.yml)
[![Test Python](https://github.com/mbruno46/chiexp/actions/workflows/python.yml/badge.svg)](https://github.com/mbruno46/chiexp/actions/workflows/python.yml)
[![Build Doc](https://github.com/mbruno46/chiexp/actions/workflows/doc.yml/badge.svg)](https://github.com/mbruno46/chiexp/actions/workflows/doc.yml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

[![Python](https://img.shields.io/badge/Python-3.10+-brightgreen.svg)](https://www.python.org)
[![Matlab](https://img.shields.io/badge/MATLAB-R2019b-brightgreen.svg)](https://www.mathworks.com/products/matlab.html)
[![Octave](https://img.shields.io/badge/Octave-5.2.0-brightgreen.svg)](https://www.mathworks.com/products/matlab.html)

A package to compute the expectation value of the chi squared defined 
from arbitrary regularizations of the inverse covariance matrix.


- **Website:** https://mbruno46.github.io/chiexp/
- **Documentation:** [web](https://mbruno46.github.io/chiexp/), [pdf][3]
- **Examples:** [python](./examples/python), [matlab](./examples/matlab)
- **Source code:** [python](./lib/python), [matlab](./lib/matlab)
- **Bug reports:** https://github.com/mbruno46/chiexp/issues


If you use this library in your publications please cite:

 - M. Bruno, R. Sommer, [On fits to correlated and auto-correlated data](https://inspirehep.net/literature/2157883) *Comput.Phys.Commun.* 285 (2023) 108643.

### Authors

Copyright (C) 2017-2022 Mattia Bruno, Rainer Sommer

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

To build the documentation locally download the git repository and 
run `make html` or `make latexpdf` from the `doc` folder. 
It requires the `sphinx-build` command so make sure that `sphinx`
is properly installed (`pip install sphinx`).
 
[1]: https://arxiv.org
[2]: https://mbruno46.github.io/chiexp
[3]: ./docs/chiexp-doc.pdf
