#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:19:16 2019

@author: tammas

This is an implementation of daisy world: a theoretical 
model of the Gaia hypothesis. This theory was initially developed by 
Lovelock (1983) to demonstrate the plausibility of living things 
interacting with, and regulating, their environment.

    Wood, A. J., G. J. Ackland, J. G. Dyke, H. T. P. Williams, and T. M. 
        Lenton, 2015: Daisywolrd: a Review. Rev. Geophys., 46, RG1001, 
        https://doi.org/10.1029/2006RG000217.
"""
if __name__=='__main__':
    print('execute main program')

import numpy as np
import matplotlib.pyplot as plt

from daisyworld import albedo
from daisyworld import temp
from daisyworld import euler
from daisyworld import beta
from daisyworld import repl

#def albedo(alphaw, alphab, alphag, aw, ab, ag):
#    """Caluclate the average planetary albedo. 
#    alphaw*aw + alphab*ab + alphag*ag
#    
#    Arguments
#    ---------
#    alphaw : float
#        Surface fraction of white daisies.
#    alphab : float
#        Surface fraction of black daisies.
#    alphag : float
#        Surface fraction of bare ground.
#    aw = 0.7 : float
#        Albedo of white daisies.
#    ab = 0.8 : float
#        Albedo of black dasies.
#    ag = 0.5 : float
#        Albedo of bare ground
#    """
#    return alphaw*aw + alphab*ab + alphag*ag
#
#
#def daisy_replicator(alpha, alphag, beta, gamma):
#    """Calculate the rate of change of a replicator (daisies). The formula is:
#    a*[ag*b-g]
#    which describes a logistic growth and constant death within a limited 
#    resource system over time t.
#    
#    Arguments
#    ---------
#    alpha :  float
#        alpha -- the fraction of surface covered by daisies.
#    alphag : float
#        alpha_g -- The fraction of available bare ground.
#    beta :  float
#        beta -- the birth rate.
#    gamma : float
#        gamma -- the death rate.
#    """
#    return alpha*(alphag*beta-gamma)
#
#
#def beta(temp, optimum=295.65, tmin=278.15, tmax=313.15, k=0.003265):
#    """Calculate the environment parameter beta. This is like a birth rate
#    that reaches a maximum at an optimal temperature, and a range of birth
#    rates is specified by a parabolic width parameter of k=0.003265.
#    1 - k*(temp-optimum)**2
#    
#    Arguments
#    ---------
#    temp : float
#        The environment temperature experienced by the daisies
#    
#    Keyword arguments
#    -----------------
#    optimum = 295.65 : float
#        The optimum temperature for daisy birthrate.
#    tmin = 278.15 : float
#        Minimum temperature for birthrate.
#    tmax = 313.15 : float
#        Maximum temperature for birthrate.
#    k = 0.003265 : float
#        Quadratic width parameter.
#    """
#    if (temp>tmin)&(temp<tmax):
#        return 1 - k*(temp-optimum)**2
#    else:
#        return 0
#
#
#def euler(initial, tendency, h=1):
#    """Integrate forward in time using Euler's method of numerical integration.
#    initial + h*tendency
#    
#    Arguments
#    ---------
#    initial :  float
#        The initial state.
#    tendency : float
#        The rate of change in the initial state.
#    
#    Keyword arguments
#    -----------------
#    h = 1 : float
#        The timestep duration.
#    """
#    return initial + h*tendency
#
#
#def local_temp(A, albedo, T, q=30):
#    """Calculate local temperature experienced by a particular daisy type or 
#    ground cover. This is a simplified version of the original.
#    q*(A-albedo)+T
#    
#    Arguments
#    ---------
#    A : float
#        Planetary albedo
#    alpha : float
#        Albedo of daisy type
#    T : float
#        Planetary temperature.
#    
#    Keyword Arguments
#    -----------------
#    q = 30 : float
#        Heat transfer coefficient
#    """
#    return q*(A-albedo)+T
#
#
#def planetary_temp(S, A, L=1.0):
#    """Calculate the planetary temperature.
#    SL(1-A) = sT**4
#    
#    Arguments
#    ---------
#    S : float
#        Incident solar energy.
#    A : float
#        Planetary albedo.
#
#    Keyword Arguments
#    -----------------
#    L = 1.0 : float
#        Normalised stellar luminosity.
#    """
#    sigma = 5.67032e-8 # Stephan-Bolzmann constant.
#    return ((S*L*(1-A))/sigma)**(1/4.)


# Define constants and variables
alphaw = 0.01  # Cover fraction of white daisies
alphab = 0.01   # Cover fraction of black daisies
p = 1           # The fraction of habitable surface
alphag = p-alphaw-alphab # Cover fraction of bare ground
aw = 0.75       # Albedo of white daisies
ab = 0.25       # Albedo of black daisies
ag = 0.5        # Albedo of bare ground
gamma = 0.3     # The death rate 
S = 1000        # Solar constant (W/m^2)
maxn = 1000     # Maximum number of iterations
tol = 0.000001  # Tollerance of solution
luminosities = np.arange(0.5,1.6, 0.002) # Stelar luminosities
alphaw_out = np.ones(len(luminosities))*np.nan # Output variable for white
alphab_out = np.ones(len(luminosities))*np.nan # Output variable for black
temp_out = np.ones(len(luminosities))*np.nan   # Output variable for temp

# Main loop for changing luminosity
for i,L in enumerate(luminosities):
    # Set a minimum for cover fractions
    if alphaw<0.01: alphaw = 0.01
    if alphab<0.01: alphab = 0.01
    alphag = p-alphaw-alphab
    # Reset counters
    n = 0
    changew, changeb = 1,1
    # Run loop for daisy earth.
    while (n<maxn) and (changew>tol) and (changeb>tol):
        # Store the initial cover fractions
        sw,sb = alphaw, alphab
        # Planetary albedo
        planet_albedo = albedo.albedo(alphaw,alphab,alphag,aw,ab,ag)
        # Planetary temperature
        T = temp.planetary_temp(S,planet_albedo, L=L)
        # Local temperature
        Tw = temp.local_temp(planet_albedo,aw,T)
        Tb = temp.local_temp(planet_albedo,ab,T)
        # Birth rate
        betaw = beta.beta(Tw)
        betab = beta.beta(Tb)
        # Change in daisies
        dawdt = repl.daisy_replicator(alphaw, alphag, betaw, gamma)
        dabdt = repl.daisy_replicator(alphab, alphag, betab, gamma)
        # Integrate
        alphaw = euler.euler(alphaw, dawdt)
        alphab = euler.euler(alphab, dabdt)
        alphag = p-alphaw-alphab
        n += 1
    # Store the output
    alphaw_out[i] = alphaw
    alphab_out[i] = alphab
    temp_out[i] = T

# Plot the results
# Cover fractions
white = plt.plot(luminosities,alphaw_out*100,'b', label='White')
black = plt.plot(luminosities,alphab_out*100,'k', label='Black')
plt.legend(loc='upper right')
plt.xlabel('Luminosity')
plt.ylabel('Surface cover %')
plt.title('Cover fractions')
plt.show()

# Planetary temperature
plt.figure()
plt.plot(luminosities,temp_out-273.15,'r')
plt.xlabel('Luminosity')
plt.ylabel('Temperature (C)')
plt.title('Planetary temperature')
plt.show()

#pyinstaller error: cannot find existing PyQt5 plugin directories
