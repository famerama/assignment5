# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:33:41 2019

@author: ru54vab
"""


for i,L in enumerate(daisyworld.luminosities):
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
        planet_albedo = albedo(alphaw,alphab,alphag,aw,ab,ag)
        # Planetary temperature
        T = temp.planetary_temp(S,planet_albedo, L=L)
        # Local temperature
        Tw = temp.local_temp(planet_albedo,aw,T)
        Tb = temp.local_temp(planet_albedo,ab,T)
        # Birth rate
        betaw = beta(Tw)
        betab = beta(Tb)
        # Change in daisies
        dawdt = daisy_replicator(alphaw, alphag, betaw, gamma)
        dabdt = daisy_replicator(alphab, alphag, betab, gamma)
        # Integrate
        alphaw = euler(alphaw, dawdt)
        alphab = euler(alphab, dabdt)
        alphag = p-alphaw-alphab
        n += 1
    # Store the output
    alphaw_out[i] = alphaw
    alphab_out[i] = alphab
    temp_out[i] = T
