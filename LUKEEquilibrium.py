# Implementation of a LUKE equilibrium.

import numpy as np
import scipy.optimize
from scipy.interpolate import RectBivariateSpline
from DREAM.Settings.LUKEMagneticField import LUKEMagneticField
from . Equilibrium import Equilibrium


class LUKEEquilibrium(Equilibrium):
    

    def __init__(self, filename):
        """
        Constructor.
        """
        self._luke_magfield = LUKEMagneticField(filename)

        self.r = self._luke_magfield.ptx[0,:]
        psi = self._luke_magfield.psi_apRp
        self.psin = (psi - psi[0]) / (psi[-1] - psi[0])

        self.R = RectBivariateSpline(
            self.psin,
            self._luke_magfield.theta,
            self._luke_magfield.ptx.T+self._luke_magfield.Rp
        )
        self.Z = RectBivariateSpline(
            self.psin,
            self._luke_magfield.theta,
            self._luke_magfield.pty.T+self._luke_magfield.Zp
        )

        self.dRdp = self.R.partial_derivative(1, 0)
        self.dRdt = self.R.partial_derivative(0, 1)
        self.dZdp = self.Z.partial_derivative(1, 0)
        self.dZdt = self.Z.partial_derivative(0, 1)

        self.lcfs = np.array([
            self._luke_magfield.ptx[:,-1] + self._luke_magfield.Rp,
            self._luke_magfield.pty[:,-1] + self._luke_magfield.Zp
        ])


    def length_to_radius(self, l, x0, nhat):
        """
        Convert the given length along a line-of-sight (originating in x0 and
        extending along nhat) to a flux surface label.

        The flux-surface label is determined by minimizing
          
            sqrt((R-R0)^2+(Z-Z0)^2),

        where (R0, Z0) is the point corresponding to the position l along the
        line-of-sight.
        """
        X = x0 + l*nhat
        R, Z = np.hypot(X[0], X[1]), X[2]

        def fun(x, R, Z):
            _R = self.R(x[0], np.mod(x[1], 2*np.pi))[0,0] - R
            _Z = self.Z(x[0], np.mod(x[1], 2*np.pi))[0,0] - Z
            return np.sqrt(_R**2 + _Z**2)


        def grad(x, R, Z):
            f = fun(x, R, Z)
            p = x[0]
            t = np.mod(x[1], 2*np.pi)
            dRdp = self.dRdp(p, t)[0,0]
            dZdp = self.dZdp(p, t)[0,0]
            dRdt = self.dRdt(p, t)[0,0]
            dZdt = self.dZdt(p, t)[0,0]

            _R = self.R(p, t)[0,0]
            _Z = self.Z(p, t)[0,0]

            return [(_R-R)*dRdp + (_Z-Z)*dZdp, (_R-R)*dRdt + (_Z-Z)*dZdt]

        # Select initial guess
        theta0 = np.arctan2(Z-self._luke_magfield.Zp, R-self._luke_magfield.Rp)
        if theta0 < 0:
            theta0 += 2*np.pi
        psi0 = (R-self._luke_magfield.Rp) / (self.R(1, theta0)[0,0] - self._luke_magfield.Rp)

        x0 = [psi0, theta0]

        # Solve
        sol = scipy.optimize.minimize(fun, x0, args=(R, Z), jac=grad, bounds=[(0, 1), (0, 2*np.pi)])

        if sol.x[0] < 0 or sol.x[0] > 1:
            #sol = scipy.optimize.root(fun, [1, 0], args=(R, Z))
            print(f'R = {R}, Z = {Z}')
            print(sol.x)
            import matplotlib.pyplot as plt
            _p = np.linspace(0, 1)
            _t = np.linspace(0, 2*np.pi)

            plt.figure()
            plt.contourf(_p, _t, np.sqrt((self.R(_p, _t)-R)**2 + (self.Z(_p, _t)-Z)**2).T)
            plt.colorbar()

            plt.figure()
            plt.plot(self.lcfs[0,:], self.lcfs[1,:], 'k', lw=2)
            plt.plot(self._luke_magfield.Rp, self._luke_magfield.Zp, 'bo')
            plt.plot(R, Z, 'rx')
            plt.show()

            raise Exception("Failed to map (R, Z) -> (psi, theta).")

        # Return corresponding DREAM radial coordinate
        return self.R(sol.x[0], 0) - self._luke_magfield.Rp


