#!/usr/bin/env python
u"""
forward.py
Written by Tyler Sutterley (05/2024)
Class for estimating non-linear parameters with a finite difference model
using a Crank-Nicolson scheme

REFERENCES:
    G. D. Egbert and S. Erofeeva, "Efficient Inverse Modeling of Barotropic
        Ocean Tides", Journal of Atmospheric and Oceanic Technology, (2002).
    F. W. Primeau and D. Newman, "Bifurcation structure of a wind-driven
        shallow water model with layer-outcropping", Ocean Modelling, (2007).

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        https://numpy.org
        https://numpy.org/docmatrix/stable/user/numpy-for-matlab-users.html
    scipy: Scientific Tools for Python
        https://docs.scipy.org/doc/
    pyproj: Python interface to PROJ library
        https://pypi.org/project/pyproj/
        https://pyproj4.github.io/pyproj/

PROGRAM DEPENDENCIES:
    arguments.py: loads nodal corrections for tidal constituents
    astro.py: computes the basic astronomical mean longitudes
    grid.py: sets up finite difference grids for tidal modeling
    utilities.py: download and management utilities for files

UPDATE HISTORY:
    Updated 05/2024: make subscriptable and allow item assignment
    Written 03/2024
"""
from __future__ import annotations

import numpy as np
import scipy.sparse.linalg
import pyTMD.solve.grid as fdgrid
import pyTMD.solve.model as model
from pyTMD.utilities import reify

class forward(model):
    """
    Class for estimating non-linear parameters using a finite
    difference model with a Crank-Nicolson scheme
    :cite:p:`Egbert:2002ge` :cite:p:`Primeau:2007ie`

    Parameters
    ----------
    grid: pyTMD.solve.grid
        Finite difference grid
    beta: float, default 0.10
        scaling factor for load tides and self attraction effects
    rho_w: float, default 1035.0
        average density of sea water [kg/m^3]
    c_drag: float, default 0.003
        non-dimensional quadratic bottom drag coefficient
    nu: float, default 100.0
        Eddy viscosity [m^2/s^2]
    Rd: float or NoneType, default None
        Rossby radius of deformation
    h2: float, default 0.609
        Love number of vertical displacement
    k2: float, default 0.302
        Love number of gravitational potential
    dt: float, default 30.0
        time step in seconds
    """
    # seconds per day
    day = 86400.0
    # # sidereal day in seconds
    # ts = 86164.0
    # 360 degrees
    turndeg = 360.0
    # degrees to radians
    deg2rad = np.pi/180.0

    # inherit model class
    def __init__(self, grid: fdgrid, **kwargs):
        # set initial model
        super().__init__(grid, **kwargs)
        # create iterator and stop index
        self.__index__ = 0
        self.__stop__ = 0

    def initialize(self, **kwargs):
        """Set the initial conditions for the model

        Parameters
        ----------
        t0: float, default 0.0
            initial time in seconds
        u0: np.ndarray, default np.zeros((ny, nx))
            initial velocity field in the x-direction
        v0: np.ndarray, default np.zeros((ny, nx))
            initial velocity field in the y-direction
        z0: np.ndarray, default np.zeros((ny, nx))
            initial height field
        k: np.ndarray, default np.zeros((ny, nx))
            friction coefficient field
        H0: np.ndarray, default np.zeros((ny, nx))
            bathymetry field
        ob: np.ndarray, default np.zeros((ny, nx), dtype=bool)
            open boundary conditions
        """
        # default initial conditions
        kwargs.setdefault('t0', 0.0)
        kwargs.setdefault('u0', np.zeros((self.ny, self.nx)))
        kwargs.setdefault('v0', np.zeros((self.ny, self.nx)))
        kwargs.setdefault('z0', np.zeros((self.ny, self.nx)))
        # initialize friction and depth fields
        kwargs.setdefault('k', np.zeros((self.ny, self.nx)))
        kwargs.setdefault('H0', np.zeros((self.ny, self.nx)))
        # open boundary conditions of the model
        kwargs.setdefault('ob', np.zeros((self.ny, self.nx), dtype=bool))
        # polynomial coefficients for calculating time-variable
        # boundary conditions
        kwargs.setdefault('pc', [0])
        # set the initial data type for all fields
        self.dtype = np.array([kwargs['u0'][0],
            kwargs['v0'][0], kwargs['z0'][0]]).dtype
        self.iscomplex = np.iscomplexobj([kwargs['u0'][0],
            kwargs['v0'][0], kwargs['z0'][0]])
        # time variable
        self.t = np.copy(kwargs['t0'])
        # model fields
        self.k = np.copy(kwargs['k'])
        self.H0 = np.copy(kwargs['H0'])
        # boundary condition polynomial coefficients
        self._pc = np.copy(kwargs['pc'])
        # allocate field for flattened open boundary condition
        self._ob = kwargs['ob'].flatten().astype(bool)
        # initial conditions for the currents and heights
        # setting data type to be uniform between fields
        u0 = kwargs['u0'].flatten().astype(self.dtype)
        v0 = kwargs['v0'].flatten().astype(self.dtype)
        z0 = kwargs['z0'].flatten().astype(self.dtype)
        # save boundary conditions
        if self.has_ob:
            self._u0 = u0[self._iob].copy()
            self._v0 = v0[self._iob].copy()
            self._z0 = z0[self._iob].copy()
        # allocate fields for the currents and heights
        if self.iscomplex:
            th = self.polynomial_sum(self._pc, self.t)
            self._u = u0.real*np.cos(th) - u0.imag*np.sin(th)
            self._v = v0.real*np.cos(th) - v0.imag*np.sin(th)
            self._z = z0.real*np.cos(th) - z0.imag*np.sin(th)
        else:
            self._u = u0.copy()
            self._v = v0.copy()
            self._z = z0.copy()
        # create links for reshaped output fields
        self.u = self._u.reshape(self.ny, self.nx)
        self.v = self._v.reshape(self.ny, self.nx)
        self.z = self._z.reshape(self.ny, self.nx)
        # state vector
        self.state = np.r_[self._u, self._v, self._z]
        # setup the model operators and the model
        self._linear_operators()
        self._factorize()

    def d0(self, M: np.ndarray):
        """Converts a vector or matrix to a sparse diagonal matrix

        Parameters
        ----------
        M: np.ndarray
            input field to be converted
        """
        return scipy.sparse.spdiags(M.flatten(), 0, M.size, M.size)

    def polynomial_sum(self, c: list | np.ndarray, t: np.ndarray):
        """
        Calculates the sum of a polynomial function of time

        Parameters
        ----------
        c: list or np.ndarray
            leading coefficient of polynomials of increasing order
        t: np.ndarray
            delta time
        """
        # convert time to array if importing a single value
        t = np.atleast_1d(t)
        return np.sum([ci * (t ** i) for i, ci in enumerate(c)], axis=0)

    def _linear_operators(self):
        """Create the operators for the finite difference model
        """
        # number of nodes
        n = self.size
        # grid spacing taking into account grid distortion
        dx, dy = self.grid.get_step_size(self.grid.x_z, self.grid.y_z)
        # identity matrix
        I = scipy.sparse.eye(n, n).tocsc()
        # identity matrix for the x and y directions
        Ix = scipy.sparse.eye(self.nx).tocsc()
        Iy = scipy.sparse.eye(self.ny).tocsc()
        # Dirichlet boundary conditions for y directions
        iin = scipy.sparse.eye(self.ny, k=1)
        iis = scipy.sparse.eye(self.ny, k=-1)
        # boundary conditions for the x directions
        if self.grid.is_global:
            # Periodic boundaries
            iix = np.arange(self.nx)
            iie = Ix[np.roll(iix, -1), :]
            iiw = Ix[np.roll(iix, 1), :]
        else:
            # Dirichlet boundaries
            iie = scipy.sparse.eye(self.nx, k=1)
            iiw = scipy.sparse.eye(self.nx, k=-1)
        # full shift identity matrices
        IE = scipy.sparse.kron(Iy, iie).tocsc()
        IW = scipy.sparse.kron(Iy, iiw).tocsc()
        IN = scipy.sparse.kron(iin, Ix).tocsc()
        IS = scipy.sparse.kron(iis, Ix).tocsc()

        # finite difference operators taking into account grid distortion
        # first order forward operators
        FDx = self.d0(1.0/dx.flatten()).dot(IE - I)
        FDy = self.d0(1.0/dy.flatten()).dot(IN - I)
        # first order backwards operators
        BDx = self.d0(1.0/dx.flatten()).dot(I - IW)
        BDy = self.d0(1.0/dy.flatten()).dot(I - IS)

        # setup flattened masks
        self.mask_z = self.grid.mask_z.flatten()
        self.mask_u = self.mask_z*IE*self.mask_z
        self.mask_v = self.mask_z*IN*self.mask_z
            
        # get the longitude and latitude of the nodes
        _, lat_z = self.grid.get_latlon(self.grid.x_z, self.grid.y_z)
        _, lat_u = self.grid.get_latlon(self.grid.x_u, self.grid.y_u)
        _, lat_v = self.grid.get_latlon(self.grid.x_v, self.grid.y_v)
        # colatitudes for nodes in radians
        th_z = self.deg2rad*(90.0 - lat_z)
        th_u = self.deg2rad*(90.0 - lat_u)
        th_v = self.deg2rad*(90.0 - lat_v)

        # Coriolis parameters at each node
        f_z = self.grid.f(th_z)
        f_u = self.grid.f(th_u)
        f_v = self.grid.f(th_v)
        # calculate Rossby deformation radius
        if (self.Rd is None):
            # default to the standard shallow water model
            # gravity at zeta nodes
            gamma_z = self.grid.gamma_0(th_z)
            # Rossby deformation radius at zeta nodes
            Rd = np.sqrt(gamma_z*self.H0)/f_z
        else:
            # Rossby deformation radius at zeta nodes
            Rd = self.Rd*np.ones((self.ny, self.nx))
        # reduced gravity at valid (wet) zeta nodes
        # taking into account load tides and self-attraction
        gp = np.zeros((self.ny, self.nx))
        ind = np.nonzero(self.H0 > 0.0)
        gp[ind] = (1.0 - self.beta) * (f_z[ind] * Rd[ind])**2 / self.H0[ind]

        # discrete approximation of the gradient
        # using backwards difference operators
        GRAD = scipy.sparse.vstack([BDx, BDy])
        # discrete approximation of the divergence terms
        # using backwards difference operators
        DIV = -GRAD.T
        DIVzu = self.d0(self.H0)*BDx*self.d0(self.mask_z)*self.d0(IE*self.mask_z)
        DIVzv = self.d0(self.H0)*BDy*self.d0(self.mask_z)*self.d0(IN*self.mask_z)

        # operators for the u nodes
        DY0 = self.d0(self.mask_z)*self.d0(IN*self.mask_z)*FDy + \
            self.d0(self.mask_z)*self.d0(1-IN*self.mask_z)*(self.d0(1.0/dy)*(-2*I)) + \
            self.d0(1-self.mask_z)*self.d0(IN*self.mask_z)*(self.d0(1.0/dy)*(2*IN))
        # 2D gradient operator (forward differences)
        GRADu = scipy.sparse.vstack([FDx, DY0])
        # 2D laplacian is the divergence of the gradient
        # using a centered difference for the laplacian
        DEL2u = DIV.dot(GRADu)
        # averaging operator for calculating the V currents at U nodes
        # which zeroes out the velocities at the masked points
        # (Dirichlet boundary conditions at the coasts)
        ISE = 0.25*(I + IE + IS + IS*IE)
        Vu = ISE*self.d0(self.mask_z)*self.d0(IN*self.mask_z)

        # operators for the v nodes
        DX0 = self.d0(self.mask_z)*self.d0(IE*self.mask_z)*FDx + \
            self.d0(self.mask_z)*self.d0(1-IE*self.mask_z)*(self.d0(1.0/dx)*(-2*I)) + \
            self.d0(1-self.mask_z)*self.d0(IE*self.mask_z)*(self.d0(1.0/dx)*(2*IE))
        # 2D gradient operator (forward differences)
        GRADv = scipy.sparse.vstack([DX0, FDy])
        # 2D laplacian is the divergence of the gradient
        # using a centered difference for the laplacian
        DEL2v = DIV.dot(GRADv)
        # averaging operator for calculating the U currents at V nodes
        # which zeroes out the velocities at the masked points
        # (Dirichlet boundary conditions at the coasts)
        INW = 0.25*(I + IN + IW + IN*IW)
        Uv = INW*self.d0(self.mask_z)*self.d0(IE*self.mask_z)

        # operator for the linear shallow water model dynamics
        # with viscous terms and coriolis effects
        # NOTE that the locations of the Coriolis terms are at
        # u nodes for the first row and v nodes for the second row
        #     __u__ __v__ __z__
        #  u |_____|_____|_____|
        #  v |_____|_____|_____|
        #  z |_____|_____|_____|
        self.M = scipy.sparse.vstack([
            scipy.sparse.hstack([-self.nu*DEL2u, -self.d0(f_u)*Vu, self.d0(gp)*FDx]),
            scipy.sparse.hstack([self.d0(f_v)*Uv, -self.nu*DEL2v, self.d0(gp)*FDy]),
            scipy.sparse.hstack([DIVzu, DIVzv, scipy.sparse.csc_matrix((n, n))])
        ]).tocsc()

    def _nonlinear_terms(self):
        """Calculate nonlinear terms for the model
        """
        # wrap data if global
        mode = 'wrap' if self.grid.is_global else 'edge'
        # translate u and v to the zeta grid
        # pad the u data and calculate the center average (at zeta nodes)
        tmp1 = np.pad(self.u, ((0, 0), (0, 1)), mode=mode)
        ubar = 0.5*(tmp1[:,:-1] + tmp1[:,1:])
        # pad the v data and calculate the center average (at zeta nodes)
        tmp2 = np.pad(self.v, ((0, 1), (0, 0)), mode='edge')
        vbar = 0.5*(tmp2[:-1,:] + tmp2[1:,:])
        # calculate the total height at the zeta nodes
        hz = self.H0 + (1.0 - self.beta) * self.z
        # calculate the (quadratic) bottom drag
        self.kappa = self.k*np.sqrt(ubar**2 + vbar**2)/hz
        tau_x = -ubar*self.kappa
        tau_y = -vbar*self.kappa
        b = np.zeros((self.shape))
        # calculate the body forces at all points
        F = np.r_[tau_x.flatten(), tau_y.flatten(), b.flatten()]
        # return the body forces at valid points
        return F[self.imask]

    def _factorize(self):
        """Initialize the linear system of equations
        """
        # number of nodes
        n = self.size
        # state vector at valid points
        self._s = self.state[self.imask]
        # identity matrix for complete state
        I = scipy.sparse.eye(3*n, 3*n).tocsc()
        # only keep rows and columns for valid points
        # using a Crank-Nicolson scheme for the time step
        A = I[self.imask,:] - (self.dt/2) * self.M[self.imask,:]
        B = I[self.imask,:] + (self.dt/2) * self.M[self.imask,:]
        # prefactor operators for inversion (LU decomposition)
        # function for solving the sparse linear system
        self.solve = scipy.sparse.linalg.factorized(A[:,self.imask])
        self.rhs = B[:,self.imask]

    def _enforce_boundary(self):
        """Enforce the model boundary conditions
        """
        # enforce open boundary conditions for currents and heights
        if self.iscomplex and self.has_ob:
            th = self.polynomial_sum(self._pc, self.t)
            self._u[self._iob] = self._u0.real*np.cos(th) - \
                self._u0.imag*np.sin(th)
            self._v[self._iob] = self._v0.real*np.cos(th) - \
                self._v0.imag*np.sin(th)
            self._z[self._iob] = self._z0.real*np.cos(th) - \
                self._z0.imag*np.sin(th)
        elif self.has_ob:
            self._u[self._iob] = self._u0.copy()
            self._v[self._iob] = self._v0.copy()
            self._z[self._iob] = self._z0.copy()

    def run(self, N: int):
        """Run the model for a number of time steps

        Parameters
        ----------
        N: int
            number of time steps to run the model
        """
        self.__stop__ = N
        try:
            self.__next__()
        except StopIteration:
            pass
        else:
            return self

    @reify
    def _iu(self) -> np.ndarray:
        """Indices of valid (wet) u nodes
        """
        return np.flatnonzero(self.mask_u)

    @reify
    def _iv(self) -> np.ndarray:
        """Indices of valid (wet) v nodes
        """
        return np.flatnonzero(self.mask_v)

    @reify
    def _iz(self) -> np.ndarray:
        """Indices of valid (wet) zeta nodes
        """
        return np.flatnonzero(self.mask_z)

    @reify
    def _iob(self) -> np.ndarray:
        """Indices of open boundary condition nodes
        """
        return np.flatnonzero(self._ob)

    @reify
    def imask(self) -> np.ndarray:
        """Indices of valid (wet) nodes in the complete state vector
        """
        return np.flatnonzero(np.r_[self.mask_u, self.mask_v, self.mask_z])

    @reify
    def iu(self) -> np.ndarray:
        """Indices of valid (wet) u nodes within the complete state vector
        """
        n = self.grid.size
        return np.flatnonzero(np.r_[self.mask_u, np.zeros((n)), np.zeros((n))])

    @reify
    def iv(self) -> np.ndarray:
        """Indices of valid (wet) v nodes within the complete state vector
        """
        n = self.grid.size
        return np.flatnonzero(np.r_[np.zeros((n)), self.mask_v, np.zeros((n))])

    @reify
    def iz(self) -> np.ndarray:
        """Indices of valid (wet) zeta nodes within the complete state vector
        """
        n = self.grid.size
        return np.flatnonzero(np.r_[np.zeros((n)), np.zeros((n)), self.mask_z])
    
    @reify
    def has_ob(self) -> np.ndarray:
        """Defines if the model has open boundary conditions
        """
        return np.any(self._ob)

    def __str__(self):
        """String representation of the ``forward`` object
        """
        properties = ['pyTMD.solve.forward']
        shape = ', '.join(map(str, self.shape))
        properties.append(f"    shape: {shape}")
        return '\n'.join(properties)

    def __iter__(self):
        """Iterate over time steps
        """
        self.__index__ = 0
        return self

    def __next__(self):
        """Calculate the fields at the next time step
        """
        # halt run if the iterator is greater than the stop index
        if (self.__index__ >= self.__stop__):
            raise StopIteration
        # calculate the non-linear (friction) terms at time step
        # to add to the right-hand side of the equation
        F = self._nonlinear_terms()
        # solve for the time step using values at previous time step
        # and the non-linear terms
        self._s = self.solve(self.rhs * self._s + self.dt * F)
        # update state vector
        self.state[self.imask] = self._s
        # add time step
        self.t += self.dt
        # update fields at time step
        self._u[self._iu] = self.state[self.iu]
        self._v[self._iv] = self.state[self.iv]
        self._z[self._iz] = self.state[self.iz]
        # enforce the boundary conditions at time step
        self._enforce_boundary()
        # add to counter
        self.__index__ += 1
        return self

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)
