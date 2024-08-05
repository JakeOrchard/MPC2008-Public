"""Simple nonlinear solvers"""

import numpy as np
import warnings
import scipy.sparse as sp


def newton_solver2(f, x0, y0=None, args=(), jac=None, tol=1E-9, maxcount=100, backtrack_c=0.5, verbose=True):
    """Simple line search solver for root x satisfying f(x)=0 using Newton direction.

    Backtracks if input invalid or improvement is not at least half the predicted improvement.

    Parameters
    ----------
    f               : function, to solve for f(x)=0, input and output are arrays of same length
    x0              : array (n), initial guess for x
    y0              : [optional] array (n), y0=f(x0), if already known
    args            : [optional] function arguments
    jac             : [optional] function that computes jacobian analytically
    tol             : [optional] scalar, solver exits successfully when |f(x)| < tol
    maxcount        : [optional] int, maximum number of Newton steps
    backtrack_c     : [optional] scalar, fraction to backtrack if step unsuccessful, i.e.
                        if we tried step from x to x+dx, now try x+backtrack_c*dx

    Returns
    ----------
    x       : array (n), (approximate) root of f(x)=0
    y       : array (n), y=f(x), satisfies |y| < tol
    """

    x, y = x0, y0
    if y is None:
        y = f(x, *args)      

    for count in range(maxcount):
        if verbose:
            printit(count, x, y)

        if np.max(np.abs(y)) < tol:
            return x, y
        
        if jac is None:
            J = obtain_J2(f, x, y, args=args)
        else:
            J = jac(x, *args)
       
        if sp.issparse(J):
            dx = sp.linalg.spsolve(J.tocsc(), -y)
        else:
            dx = np.linalg.solve(J, -y)

        # backtrack at most 29 times
        for bcount in range(30):
            try:
                ynew = f(x + dx, *args)
            except ValueError:
                if verbose:
                    print('backtracking\n')
                dx *= backtrack_c
            else:
                predicted_improvement = -np.sum((J @ dx) * y) * ((1 - 1 / 2 ** bcount) + 1) / 2
                actual_improvement = (np.sum(y ** 2) - np.sum(ynew ** 2)) / 2
                if actual_improvement < predicted_improvement / 2:
                    if verbose:
                        print('backtracking\n')
                    dx *= backtrack_c
                else:
                    y = ynew
                    x += dx
                    break
        else:
            raise ValueError('Too many backtracks, maybe bad initial guess?')
    else:
        raise ValueError(f'No convergence after {maxcount} iterations')


def obtain_J2(f, x, y, args=(), h=1E-5):
    """Finds Jacobian f'(x) around y=f(x)"""
    nx = x.shape[0]
    ny = y.shape[0]
    J = np.empty((ny, nx))

    for i in range(nx):
        dx = h * (np.arange(nx) == i)
        J[:, i] = (f(x + dx, *args) - y) / h
    return J