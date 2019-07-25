# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:44:06 2019

@author: ru54vab
"""

def euler(initial, tendency, h=1):
    """Integrate forward in time using Euler's method of numerical integration.
    initial + h*tendency
    
    Arguments
    ---------
    initial :  float
        The initial state.
    tendency : float
        The rate of change in the initial state.
    
    Keyword arguments
    -----------------
    h = 1 : float
        The timestep duration.
    """
    return initial + h*tendency