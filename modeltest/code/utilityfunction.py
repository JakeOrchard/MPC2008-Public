#%%
import numpy as np
from numba import njit
import scipy.sparse as sp

@njit
def f_invuc_dc(uc, dc, xi=1, eis=1, psi=0.5):
    
    if xi==1:
        c = (psi / uc) ** eis * dc ** ((1 - psi) * (eis - 1))

    else:
        cbase =     psi   ** (1/xi) * 1 ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * dc ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        c  = (psi)**(eis/xi) * uc ** (-eis) *  complement ** ((eis - xi) / xi)

    return c

@njit
def f_uc_dc(c, dc, xi=1, eis=1, psi=0.5):
    
    if xi==1:
        marguc =  psi * c ** (- 1 / eis) * dc ** ((1 - psi) * (1 - 1/eis))

    else:
        cbase =     psi   ** (1/xi) * 1 ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * dc ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        marguc = psi ** (1 / xi) * c **(- 1 / eis) * complement ** ((1 - xi / eis) / xi)

    return marguc  

@njit
def f_u(c, d, xi=1, eis=1, psi=0.5):
    
    if eis==1 and xi==1:
        u =  psi * np.log(c) + (1 - psi) * np.log(d)
    elif eis!=1:
        u =  (c ** (psi * (1 - 1/eis)) * d ** ((1 - psi) * (1 - 1/eis)) - 1) / (1 - 1/eis)
    else:
        cbase =     psi   ** (1/xi) * c ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * d ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        u = (complement ** (1 - 1 / eis) - 1) / (1 - 1 / eis)

    return u       

if __name__=='__main__':
    for xi in [0.5,1,2]:
        for eis in [0.5,1,2]:
            for psi in [0.25,0.5,0.75]:
                for dc in [0.5,1,2]:
                    for c in [0.1,0.2,1]:
                        f_uc(1, 1, xi=xi, eis=eis, psi=psi)



def f_uc(c, d, xi=1, eis=1, psi=0.5):
    
    if xi==1:
        marguc =  psi * c ** (psi * (1 - 1/eis) - 1) * d ** ((1 - psi) * (1 - 1/eis))

    else:
        cbase =     psi   ** (1/xi) * c ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * d ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        marguc = (psi / c) ** (1 / xi) * complement ** ((1 - xi / eis) / xi)

    return marguc  

def f_ud(c, d, xi=1, eis=1, psi=0.5):
    
    if xi==1:
        margud =  (1-psi) * c ** (psi * (1 - 1/eis)) * d ** ((1 - psi) * (1 - 1/eis) - 1)

    else:
        cbase =     psi   ** (1/xi) * c ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * d ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        margud = ((1-psi) / d) ** (1 / xi) * complement ** ((1 - xi / eis) / xi)

    return margud       


def f_invuc(marguc, d, xi=1, eis=1, psi=0.5):
    
    if xi==eis:
        c =  psi / marguc
        c = (psi ** (1 / xi) / marguc) ** eis 

    elif xi==1:
        c =  (marguc / psi) ** (1 / (psi * (1 - 1/eis) - 1)) * d ** ( - (1 - psi) * (1 - 1/eis) / (psi * (1 - 1/eis) - 1))

    else:        
        raise NotImplementedError
        # newton solve
        # f_uc(c, d, xi, eis, psi) - marguc
        #     marguc - (psi / c) ** (1 / xi) * complement ** ((1 - xi / eis) / xi)

    return c 


def f_invucc(marguc, d, xi=1, eis=1, psi=0.5):
    
    if xi!=1:
        cbase =     psi   ** (1/xi) * c ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * d ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** (xi / (xi -1))

        cshare = cbase / complement_base

        # verify sign
        jac = ( - 1 / xi * psi ** (1 / xi) * c ** (-1 / xi - 1) * 
                (1 - xi * (1 - xi / eis) / (xi -1) * cshare) 
                * complement )

    jac = sp.diags(jac.reshape([np.prod(np.shape(marguc)),]))

    return jac           


def f_ucc(c, uc, d, xi=1, eis=1, psi=0.5):
    c = c.reshape(np.shape(uc))
    
    if xi==1:
        marguc = f_uc(c, uc, d, xi=xi, eis=eis, psi=psi)

        jac =  (psi * (1 - 1/eis) - 1) * marguc / c

    else:
        cbase =     psi   ** (1/xi) * c ** ((xi -1) / xi)
        dbase = (1 - psi) ** (1/xi) * d ** ((xi -1) / xi)
        complement_base = cbase + dbase
        complement =  complement_base ** ((1 - xi / eis) / (xi -1))

        cshare = cbase / complement_base

        jac = ( - 1 / xi * psi ** (1 / xi) * c ** (-1 / xi - 1) * 
                (1 - xi * (1 - xi / eis) / (xi -1) * cshare) 
                * complement )

    jac = sp.diags(jac.reshape([np.prod(np.shape(c)),]))


    return jac        
          

# def f_resid(c, uc, d, xi, eis, psi):
#     c = c.reshape(np.shape(uc))

#     resid = f_uc(c, uc, d, xi, eis, psi) - uc

#     resid = resid.reshape([np.prod(np.shape(c)),])

#     return resid    
