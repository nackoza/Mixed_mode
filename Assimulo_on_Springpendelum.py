# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 09:46:49 2021

@author: x-mikael.gustafsson
"""
import numpy as np
from assimulo.solvers import ImplicitEuler, ExplicitEuler, RungeKutta4
from assimulo.problem import Explicit_Problem
from Springpendelum import f


y_0 = [6,0,0,np.pi/4] 
t_0 = 0
tfinal = 5
prob = Explicit_Problem(f,y_0,t_0)

sim = ImplicitEuler(prob)
sim.h = 0.001
h = sim.h
t, y = sim.simulate(tfinal)
sim.plot()

    
I = [t_0, tfinal]








































