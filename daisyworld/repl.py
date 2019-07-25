# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:41:22 2019

@author: ru54vab
"""

def daisy_replicator(alpha, alphag, beta, gamma):
    """Calculate the rate of change of a replicator (daisies). The formula is:
    a*[ag*b-g]
    which describes a logistic growth and constant death within a limited 
    resource system over time t.
    
    Arguments
    ---------
    alpha :  float
        alpha -- the fraction of surface covered by daisies.
    alphag : float
        alpha_g -- The fraction of available bare ground.
    beta :  float
        beta -- the birth rate.
    gamma : float
        gamma -- the death rate.
    """
    return alpha*(alphag*beta-gamma)