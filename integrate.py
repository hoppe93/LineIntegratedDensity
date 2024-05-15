
import matplotlib.pyplot as plt
import numpy as np
import warnings
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator


def integrate(ne, r, equilibrium, x0, nhat, epsrel=1e-4, l1=None, l2=None):
    """
    Calculate the line-integrated density nBar corresponding to the given
    electron density ``ne``, given as a function of radius ``r``. The
    calculation is to be carried out within the given magnetic equilibrium.
    Integration commences from the given origin ``x0`` in the direction
    ``nhat``.

    :param ne:          Electron density to integrate.
    :param r:           Radius vector corresponding to ne.
    :param equilibrium: Magnetic equilibrium to use for the integration.
    :param x0:          Origin of diagnostic line-of-sight.
    :param nhat:        Direction of diagnostic line-of-sight.
    :param epsrel:      Relative tolerance for integration.
    :param l1:          First point along line-of-sight which intersects the plasma. Automatically calculated if not specified.
    :param l2:          Second point along line-of-sight which intersects the plasma. Automatically calculated if not specified.
    """
    x0 = np.array(x0)
    nhat = np.array(nhat)

    nefunc = PchipInterpolator(x=r, y=ne)

    # Calculate integration limits
    if l1 is None or l2 is None:
        l1, l2 = equilibrium.find_intersections(x0, nhat)

    if l1 is None or l2 is None:
        warnings.warn(f'No intersection of line-of-sight with plasma found.')
        return 0

    _l1 = min(l1, l2)
    _l2 = max(l1, l2)
    l1, l2 = _l1, _l2

    # Perform integration
    nBar, _ = quad(lambda l : nefunc(equilibrium.length_to_radius(l, x0, nhat)), l1, l2, limit=200, epsabs=0, epsrel=epsrel)

    return nBar


