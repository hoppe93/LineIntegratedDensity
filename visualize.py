
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator


def visualize(equilibrium, x0, nhat):
    """
    Visualize the line-of-sight.
    """
    # Plot LCFS
    theta = np.linspace(0, 2*np.pi, 200)
    fig, ax = plt.subplots(1, 1, figsize=(6, 8))

    R, Z = equilibrium.lcfs[0,:], equilibrium.lcfs[1,:]
    ax.plot(R, Z, 'k', lw=2)
    ax.axis('equal')

    # Plot line-of-sight
    l = np.linspace(0, 2*max(Z))
    l1, l2 = equilibrium.find_intersections(x0, nhat)

    lX = lambda l : np.sqrt((x0[0] + l*nhat[0])**2 + (x0[1] + l*nhat[1])**2)
    lZ = lambda l : x0[2] + l*nhat[2]

    ax.plot(lX(l), lZ(l), 'tab:blue', lw=2)
    if l1 is not None:
        print(f'Intersection at l1 = {l1}.')
        ax.plot(lX(l1), lZ(l1), 'rx', ms=8, mew=3)
    else:
        print('No l1 intersection point.')
    if l2 is not None:
        print(f'Intersection at l2 = {l2}.')
        ax.plot(lX(l2), lZ(l2), 'rx', ms=8, mew=3)
    else:
        print('No l2 intersection point.')

    fig.tight_layout()


def plotDensity(ne, r, equilibrium, x0, nhat):
    """
    Plot the electron density as a function of the line-of-sight
    distance l.
    """
    x0 = np.array(x0)
    nhat = np.array(nhat)

    nefunc = PchipInterpolator(x=r, y=ne)
        
    # Evaluate intersection points
    l1, l2 = equilibrium.find_intersections(x0, nhat)

    l = np.linspace(l1, l2, 400)
    _r  = np.zeros((l.size,))
    _ne = np.zeros((l.size,))
    for i in range(l.size):
        _r[i]  = equilibrium.length_to_radius(l[i], x0, nhat)
        _ne[i] = nefunc(equilibrium.length_to_radius(l[i], x0, nhat))

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    axs[0].plot(r, ne)
    axs[0].set_xlabel('$r$')
    axs[0].set_ylabel('$n_e$')

    axs[1].plot(l, _r)
    axs[1].set_xlabel('$l$')
    axs[1].set_ylabel('$r$')

    axs[2].plot(l, _ne)
    axs[2].set_xlabel('$l$')
    axs[2].set_ylabel('$n_e$ (m$^{-3}$)')

    fig.tight_layout()


