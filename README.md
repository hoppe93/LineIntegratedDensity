# Synthetic line-integrated density diagnostic for DREAM
This repository cotnains a synthetic diagnostic for line-integrated density
measurements, designed specifically for
[DREAM](https://github.com/chalmersplasmatheory/DREAM).

## Example use
The synthetic diagnostic can be used in the following way:

```python3
from DREAM import DREAMOutput
import LineIntegratedDensity as LID
import matplotlib.pyplot as plt
import numpy as np
import sys


def main():
    do = DREAMOutput('output.h5')
    eq = LID.LUKEEquilibrium(do.settings.radialgrid.num_filename)

    # Diagnostic position
    x0 = [.9030, 0, 0.65]
    # Diagnostic line-of-sight
    nhat = [0, 0, -1]

    r = do.grid.r[:]
    t = do.grid.t[:]
    nel = np.zeros((t.size,))
    l1, l2 = eq.find_intersections(x0, nhat)

    for i in range(t.size):
        ne = do.eqsys.n_cold[i,:]
        nel[i] = LID.integrate(ne, r, eq, x0, nhat, l1=l1, l2=l2)

    plt.plot(t, nel)
    plt.xlabel('Time')
    plt.ylabel('Line-integrated density (m$^{-2}$)')
    plt.show()


if __name__ == '__main__':
    sys.exit(main())
```

