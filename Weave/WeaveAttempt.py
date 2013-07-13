# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 22:59:12 2013

@author: David

Attempt to see if weave.inline works
"""

from scipy.weave import converters
from scipy import weave
import numpy

class Grid:
    """A simple grid class that stores the details and solution of the
    computational grid."""
    def __init__(self, nx=10, ny=10, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.u = numpy.zeros((nx, ny), 'd')
        # used to compute the change in solution in some of the methods.
        self.old_u = self.u.copy()

    def setBCFunc(self, func):
        """Sets the BC given a function of two variables."""
        xmin, ymin = self.xmin, self.ymin
        xmax, ymax = self.xmax, self.ymax
        x = numpy.arange(xmin, xmax + self.dx*0.5, self.dx)
        y = numpy.arange(ymin, ymax + self.dy*0.5, self.dy)
        self.u[0 ,:] = func(xmin,y)
        self.u[-1,:] = func(xmax,y)
        self.u[:, 0] = func(x,ymin)
        self.u[:,-1] = func(x,ymax)

    def computeError(self):
        """Computes absolute error using an L2 norm for the solution.
        This requires that self.u and self.old_u must be appropriately
        setup."""
        v = (self.u - self.old_u).flat
        return numpy.sqrt(numpy.dot(v,v))


class LaplaceSolver:
    """A simple Laplacian solver that can use different schemes to
    solve the problem."""
    def __init__(self, grid, stepper='numeric'):
        self.grid = grid
        self.setTimeStepper(stepper)

    def slowTimeStep(self, dt=0.0):
        """Takes a time step using straight forward Python loops."""
        g = self.grid
        nx, ny = g.u.shape
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u

        err = 0.0
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                tmp = u[i,j]
                u[i,j] = ((u[i-1, j] + u[i+1, j])*dy2 +
                          (u[i, j-1] + u[i, j+1])*dx2)*dnr_inv
                diff = u[i,j] - tmp
                err += diff*diff

        return numpy.sqrt(err)
    
    def inlineTimeStep(self, dt=0.0):
        """Takes a time step using inlined C code -- this version uses
        blitz arrays."""
        g = self.grid
        nx, ny = g.u.shape
        dx2, dy2 = g.dx**2, g.dy**2
        dnr_inv = 0.5/(dx2 + dy2)
        u = g.u

        code = """
               #line 120 "laplace.py" (This is only useful for debugging)
               double tmp, err, diff;
               err = 0.0;
               for (int i=1; i<nx-1; ++i) {
                   for (int j=1; j<ny-1; ++j) {
                       tmp = u(i,j);
                       u(i,j) = ((u(i-1,j) + u(i+1,j))*dy2 +
                                 (u(i,j-1) + u(i,j+1))*dx2)*dnr_inv;
                       diff = u(i,j) - tmp;
                       err += diff*diff;
                   }
               }
               return_val = sqrt(err);
               """
        # compiler keyword only needed on windows with MSVC installed
        err = weave.inline(code,
                           ['u', 'dx2', 'dy2', 'dnr_inv', 'nx', 'ny'],
                           type_converters=converters.blitz,
                           compiler = 'gcc')
        return err

    def setTimeStepper(self, stepper='numeric'):
        """Sets the time step scheme to be used while solving given a
        string which should be one of ['slow', 'numeric', 'blitz',
        'inline', 'fastinline', 'fortran']."""
        if stepper == 'slow':
            self.timeStep = self.slowTimeStep
        elif stepper == 'inline':
            self.timeStep = self.inlineTimeStep
        else:
            self.timeStep = self.numericTimeStep

    def solve(self, n_iter=0, eps=1.0e-16):
        err = self.timeStep()
        count = 1

        while err > eps:
            if n_iter and count >= n_iter:
                return err
            err = self.timeStep()
            count = count + 1

        return count