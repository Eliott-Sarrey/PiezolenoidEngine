import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.misc import derivative
#import tikzplotlib
import DataFile.py

om = np.sqrt(4*k*m - c**2)/(2*m)
T = 2*pi/om
print(T)


Am = pi*Rm**2
Ag = pi*Rg**2

c = 8*pi*eta*L*Am**2/Ag**2
R = pi*2*Rs*N


def I(x, t, dxdt):
    
    #print(dphidt(x, t, dxdt))
      
    if 0 < t%T < d*T/2:
        return(Imax + dphidt(x, t, dxdt))
    if T/2 < t%T < (1+d)*T/2:
        return(-Imax + dphidt(x, t, dxdt))
    else:
        return(dphidt(x, t, dxdt))
    
def V(x, t, dxdt):
    return( I(x, t, dxdt) / R)
    
def I2(t):
    return (4*Imax / pi) * np.sum([1 / (2*n+1) * np.cos((2*n+1)*t*2*pi/T)*np.sin(n*pi + (2*n+1)*pi*d/2)*np.cos(n*pi) for n in range(1, 10)])

def F(x, t, dxdt):
    z = x+5e-2
    return( pi*Rm**2*M0*mu0*N*I(z, t, dxdt) / (2*Ls) * ( (z  + Lm) / np.sqrt((z + Lm)**2 + Rs**2) + (z - Ls) / np.sqrt((z - Ls)**2 + Rs**2) - (z ) / np.sqrt((z )**2 + Rs**2) - (z + Lm - Ls) / np.sqrt((z + Lm - Ls)**2 + Rs**2)) )
    
def dphidt(x, t, dxdt):
    z = x+5e-2
    return( dxdt*pi*Rs**2*M0*mu0*N / (2*Lm) * ( (z + Ls) / np.sqrt((z + Ls)**2 + Rm**2) +  (z - Lm) / np.sqrt((z - Lm)**2 + Rm**2) -  z / np.sqrt((z )**2 + Rm**2) - (z + Ls - Lm) / np.sqrt((z + Ls - Lm)**2 + Rm**2)) )
    

def ode_sys(t, X):
        x=X[0]
        dx_dt=X[1]
        d2x_dt2=-c*dx_dt/m - k*x / m + F(x, t, dx_dt)
        return [dx_dt, d2x_dt2]