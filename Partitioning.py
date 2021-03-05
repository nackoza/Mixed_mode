# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 16:02:01 2021

@author: x-mikael.gustafsson
"""
import numpy as np

from Springpendelum import f 
from Assimulo_on_Springpendelum import y ,h 



class Partitioning:
    
    def __init__( self , y , h ):
        
        self.y = y    ### solution points 
        self.h = h
        self.n = y.shape[ 0 ]
        self.m = y.shape[ 1 ]
        
        
    def Make_R( self , E , A ):
        N = A.shape[0]
        I = np.eye( N )
        P = np.zeros( ( N , N ) )
        for i in range( N ):
            if E[ i ] == 1:
                P[ i , i ] = 1
            else:
                P[ i , i ] = 0
        B =  I - self.h * np.dot( ( I - P ) , A ) 
        R = np.dot( np.linalg.inv( B ) , ( I + self.h * np.dot( P , A ) ) )
        return R
    
    def Make_E( self ):
        E = [ 1 for i in range( self.m ) ]
        Set  = [ i for i in range( self.m ) ]
        for i in range( self.n ) :
            A = self.Jacobian( y[ i , : ] )    
            S = np.linalg.eigvals( A ).real
            
            if all( S < 0 ):        
                R = self.Make_R( E , A )
                S_1 = np.linalg.eigvals( R )
                
                if all( abs( S_1 ) <= 1 ):
                    pass
                else:
                    for k in Set:
                        if abs( S_1[ k ] ) > 1:                    
                            E[ k ] = 0
                            Set.remove( k )
                            
        return E
    
    def Jacobian(self,x):
        eps = 1e-5
        n = np.shape( x )[ 0 ]
        A = np.zeros( ( n , n ) )
        h = np.zeros( n )
        for i in range( n ):
            h[ i ] = eps
            A[ i , : ] = ( f( 1 , x + h ) - f( 1 , x ) ) / eps
            h[ i ] = 0
        return A

parti = Partitioning(y,h)

E = parti.Make_E()  # Whant to choose Explicit and implicit method

print(E,'0 means "fast" ')












































































