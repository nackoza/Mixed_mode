# -*- cxoding: utf-8 -*-
"""
Created on Mon Feb  8 13:49:29 2021

@author: x-mikael.gustafsson
"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def f(t,x):
    x_1 = x[0]
    x_2 = x[1]
    x_3 = x[2]
    x_4 = x[3]
    l_0 = 1
    g = 9.81
    m = 3
    k = 10
    A = np.array([x_2, (l_0 + x_1)*(x_4**2) - (k/m)*x_1 + g*np.cos(x_3) ,x_4 , -(g /(l_0 + x_1))*np.sin(x_3) - (2*x_2 /(l_0 + x_1))*x_4])
    return A

























