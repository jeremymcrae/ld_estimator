## ld_estimator: A package to estimate linkage disequilibrium
Computes linkage disequilibrium.

### Installation
The simplest way to install ld_estimator is through pip:
```sh
pip install ld_estimator
```

### Usage
Use ld_estimator within a python environment
```python
from ld_estimator import pairwise_ld

var1 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (0, 0), (1, 1)]
var2 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (1, 1), (1, 1)]
is_haploid = [False, False, False, False, False, False, False, True, True, True]
ld = pairwise_ld(var1, var2, is_haploid)
print(ld.dprime)
print(ld.r_squared)
```
