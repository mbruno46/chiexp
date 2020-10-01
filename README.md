# ChiExp

A package to compute the expectation value of the chi squared defined 
from arbitrary regularizations of the inverse covariance matrix, following the results of 
Ref. [[1][1]].


## Getting started

The package is written **primarily** for Matlab 
and can be imported by typing

```matlab
>> addpath '/path/to/chiexp/directory/lib/matlab'
>> help chiexp
```

but we provide also a Python package contained
in the directory `lib/python`, that can be easily 
imported

```python
>>> import sys
>>> sys.path.append('/path/to/chiexp/directory/lib/python/chiexp')
>>> from chiexp import *
>>> help(chisquare)
```

## Documentation

Examples for Matlab can be found in the folder [example](./examples/)
and a full documentation can be found in [HTML][2]
or in [PDF][3].

If you use this library in your publications please consider citing:

 - [[1][1]] M. Bruno, R. Sommer, In preparation.

[1]: https://arxiv.org
[2]: https://mabruno.gitlab.io/chiexp
[3]: ./doc/chiexp-doc.pdf
