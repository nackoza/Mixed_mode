# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 08:43:34 2021
@author: x-mikael.gustafsson
"""
import numpy as np
from scipy import optimize
from Springpendelum import f  
import matplotlib.pyplot as plt
from Partitioning import E , h
from Assimulo_on_Springpendelum import  y , y_0 , I



class Mixed_Mode:
    
    def __init__( self , E , I , h ):
        self.E = E  ## This gives which variables must be considered as "fast"
        self.lenE = len(E)
        self.I = I  ## This is the interval 
        self.h = h  ## This is the step size, must be the same as we partitioned for
        self.N = ( I[ 1 ] - I[ 0 ]  ) / self.h 

    def simulate( self , y_0 ):
        y_n = y_0.copy()
        solution = y_0
        for i in range( int( self.N ) + 1 ):
            y_np1 = self.mixed_mode( y_n )
            solution = np.vstack( ( solution , y_np1 ) )
            y_n = y_np1.copy()
        
        return solution
            
    def mixed_mode( self , y_n ):
        "Firstly we split variables into slow and fast, I try to make it readable"
        y_S_n = np.array( [ y_n[ i ] for i in range(self.lenE) if self.E[i] == 1 ] )
        y_F_n = np.array( [ y_n[ i ] for i in range(self.lenE) if self.E[i] == 0 ] )
        "We apply a Explicit solver on the slow variables"
        y_S_np1 = self.Explicit_step( y_S_n , y_n )
        "We apply a implicit solver on the fast variables"
        y_F_np1 = self.Implicit_step( y_S_np1 , y_F_n )
        "We collect the variables into a vector y_np1"
        y_np1 = self.Collect( y_S_np1 , y_F_np1 )
        
        return y_np1
    
    def Collect( self , y_S_np1 , y_F_np1 ):
        w = []
        k = 0
        j = 0 
        for i in range( len( self.E ) ):
            if self.E[i] == 1:
                w.append( y_S_np1[ k ] )
                k = k + 1
            else:
                w.append( y_F_np1[ j ] )
                j = j + 1
        
        return np.array(w)
    
    def Explicit_step( self , y_S_n , y_n ):
        y_S_np1 = y_S_n + self.h * self.f_S(y_n)
        
        return y_S_np1

    def Implicit_step( self , y_S_np1 , y_F_n ):       
        def g( y_F ):
            f_F = lambda y_F: f( 1 , self.Collect( y_S_np1 , y_F ) )
            yf = np.array( [ f_F(y_F)[ i ] for i in range( len( self.E ) ) if self.E[ i ] == 0 ] ) 
            
            return - y_F + y_F_n + self.h * yf
        y_F_np1 = optimize.fsolve( g , y_F_n )
        
        return y_F_np1
        
    def f_S( self , y_n ):
        "We need to evaluate the whole f and then get the slow values, NOT OPTIMAL"
        f_i = f( 1 , y_n ) 
        ys = np.array( [ f_i[ i ] for i in range( len( self.E ) ) if self.E[ i ] == 1 ] )
        
        return ys

 

sim = Mixed_Mode(E,I,h)

sol = sim.simulate(y_0)

N = int(5 / h)
x_axis = np.linspace(0,5,N+1)


error1 = []
error2 = []
error3 = []
error4 = []
for j in range(N+1):
    
    error1.append(np.linalg.norm( sol[j,0] - y[j,0]) )
    error2.append(np.linalg.norm( sol[j,1] - y[j,1]) )
    error3.append(np.linalg.norm( sol[j,2] - y[j,2]) )
    error4.append(np.linalg.norm( sol[j,3] - y[j,3]) )
    

fig, ax = plt.subplots()
#ax.plot(x_axis,sol)
line1 = ax.plot(x_axis,error1,label="test1")
line2 = ax.plot(x_axis,error2,label="test2")
line3 = ax.plot(x_axis,error3,label="test3")
line4 = ax.plot(x_axis,error4,label="test4")

plt.show()


































