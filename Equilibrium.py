# Helper class for equilibrium

import numpy as np


class Equilibrium:
    

    def __init__(self):
        """
        Constructor.
        """
        pass


    def find_intersections(self, X0, nhat):
        """
        Find the parameters l1 and l2 at which the line-of-sight originating at
        x0, and directed along nhat, intersects the last closed flux surface of
        this equilibrium.
        """
        l1, l2 = None, None
        for i in range(self.lcfs.shape[1]):
            x0, z0 = self.lcfs[:,i]

            if i == self.lcfs.shape[1]-1:
                x1, z1 = self.lcfs[:,0]
            else:
                x1, z1 = self.lcfs[:,i+1]

            _l1, _l2 = self._los_intersection(X0, nhat, (x0, z0), (x1, z1))
            if _l1 is not None: l1 = _l1
            if _l2 is not None: l2 = _l2

            if l1 is not None and l2 is not None:
                break

        return l1, l2


    def _los_intersection(self, X0, nhat, p0, p1):
        """
        Check if the line-of-sight originating in x0, directed along nhat,
        intersects the line extending between p0=(x0,z0) and p1=(x1,z1).
        """
        xi,  zi  = p0[0], p0[1]
        xi1, zi1 = p1[0], p1[1]

        x0, y0, z0 = X0[0], X0[1], X0[2]
        nx, ny, nz = nhat[0], nhat[1], nhat[2]

        a0 = xi**2 - x0**2 - y0**2 + 2*xi*(z0-zi)*(xi1-xi)/(zi1-zi) +\
             (z0-zi)**2 * ((xi1-xi)/(zi1-zi))**2

        a1 = 2*xi*nz*(xi1-xi)/(zi1-zi) + 2*nz*(z0-zi)*((xi1-xi)/(zi1-zi))**2 -\
             2*(x0*nx + y0*ny)

        a2 = nz**2 * ((xi1-xi)/(zi1-zi))**2 - nx**2 - ny**2

        # Check for real solution
        sqr = -a0/a2 + a1**2/(4*a2**2)
        if sqr < 0:
            return None, None

        # Obtain the two roots
        l1 = -a1/(2*a2) + np.sqrt(sqr)
        l2 = -a1/(2*a2) - np.sqrt(sqr)

        # Check for 0 < t < 1
        t1 = (z0-zi+l1*nz)/(zi1-zi)
        t2 = (z0-zi+l2*nz)/(zi1-zi)

        if t1 < 0 or t1 > 1:
            l1 = None
        if t2 < 0 or t2 > 1:
            l2 = None

        return l1, l2


    def length_to_radius(self, l, x0, nhat):
        """
        Convert the given length along a line-of-sight (originating in x0 and
        extending along nhat) to a flux surface label.
        """
        raise NotImplementedError("The method 'length_to_radius()' has not been implemented for this equilibrium type.")


