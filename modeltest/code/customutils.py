from numba import njit

@njit
def setmin3(x, xmin):
    """Set 3-dimensional array x where each row is ascending equal to equal to max(x, xmin)."""
    ni, nj, nk = x.shape
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                if x[i, j, k] < xmin:
                    x[i, j, k] = xmin
                else:
                    break