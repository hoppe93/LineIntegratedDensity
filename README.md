# Synthetic line-integrated density diagnostic for DREAM
This repository cotnains a synthetic diagnostic for line-integrated density
measurements, designed specifically for
[DREAM](https://github.com/chalmersplasmatheory/DREAM).

## Example use
The synthetic diagnostic has a number of interfaces which are suitable in
different cases.

### integrate()
Calculate the line-integrated density corresponding to a given density profile.
Input parameters:
- ``ne``: Vector of densities to integrate.
- ``r``: Radial grid corresponding to ``ne``.
- ``equilibrium``: Equilibrium object to do the integration in.
- ``x0``: Detector position.
- ``nhat``: Detector viewing direction.
- ``epsrel``: Relative tolerance to use for integration.
- ``l1``: First point of intersection between line-of-sight and plasma. Optional; if not provided it is calculated automatically.
- ``l2``: Second point of intersection between line-of-sight and plasma. Optional; if not provided it is calculated automatically.

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

### integrate_dream()
Calculate the line-integrated density for a given subset of (or all) time
points in a ``DREAMOutput`` object.

Input parameters:
- ``do``: The ``DREAMOutput`` object to integrate.
- ``x0``: Detector position.
- ``nhat``: Detector viewing direction.
- ``nthreads``: Number of threads to parallelize over.
- ``t``: Subset of time indices to integrate.

```python3
from DREAM import DREAMOutput
import LineIntegratedDensity as LID
import matplotlib.pyplot as plt
import numpy as np
import sys


def main():
    do = DREAMOutput('output.h5')

    # Diagnostic position
    x0 = [.9030, 0, 0.65]
    # Diagnostic line-of-sight
    nhat = [0, 0, -1]

    nel = LID.integrate_dream(do, x0, nhat, nthreads=8)

    plt.plot(do.grid.t, nel)
    plt.xlabel('Time')
    plt.ylabel('Line-integrated density (m$^{-2}$)')
    plt.show()


if __name__ == '__main__':
    sys.exit(main())
```

