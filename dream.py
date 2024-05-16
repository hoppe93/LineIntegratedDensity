# DREAM API

import numpy as np
from .integrate import integrate
from .LUKEEquilibrium import LUKEEquilibrium
import DREAM.Settings.RadialGrid as RadialGrid

from multiprocessing import Pool


def integrate_dream(do, x0, nhat, nthreads=1, t=None):
    """
    Evaluate the line-integrated density for the specified subset of time steps
    in the given DREAMOutput object.
    """
    if do.settings.radialgrid.type == RadialGrid.TYPE_NUMERICAL:
        eq = LUKEEquilibrium(do.settings.radialgrid.num_filename)
    else:
        raise Exception(f"Cannot operate in an equilibrium of type '{do.settings.radialgrid.type}'. Only numerical equilibria are supported.")

    if t is None:
        tv = do.grid.t[:]
        t = range(do.grid.t.size)
    else: tv = do.grid.t[t]

    r = do.grid.r[:]
    ne = do.eqsys.n_cold[:]
    nel = np.zeros((tv.size,))
    l1, l2 = eq.find_intersections(x0, nhat)
    nea = [ne[i,:] for i in t]

    N = len(nea)
    intg = Integrand(r=r, eq=eq, x0=x0, nhat=nhat, l1=l1, l2=l2)

    if nthreads == 1:
        for i in range(len(nea)):
            nel[i] = intg(nea[i])
            print(f'Time step {i+1}/{N}')
    else:
        with Pool(nthreads) as p:
            for i, lid in enumerate(p.imap(intg, nea)):
                nel[i] = lid
                print(f'\rTime step {i+1}/{N}', end='')

        print('\n', end='')

    return nel


class Integrand:
    

    def __init__(self, r, eq, x0, nhat, l1, l2):
        self.r = r
        self.eq = eq
        self.x0 = x0
        self.nhat = nhat
        self.l1 = l1
        self.l2 = l2


    def __call__(self, ne):
        return integrate(
            ne, self.r, equilibrium=self.eq,
            x0=self.x0, nhat=self.nhat,
            l1=self.l1, l2=self.l2
        )


