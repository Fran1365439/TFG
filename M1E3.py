# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 16:10:09 2019

@author: Fran
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def F(x,y,z):
    return z

def G(x,y,z):
    return y+np.e**(x)


def solucionExacta(x):
    """
    Solución exácta del problema evaluada en el punto x.
    El valor de s se ha de definir dentro.
    """
    s = 0.96842369
    return (np.e**x/4)*(2*x-2*s-1)+(1/4)*np.e**(2*s-x)

def P1s(s, imprimir):
    """
    Resolvemos el problema P1s con el método de RK3 de derecha a izquierda.
    Denotamos a la solución por u
    """
    a, b = s, 0
    N = 20
    h = (b-a)/N
    
    x = np.zeros(21, dtype = np.double)
    y = np.zeros(21, dtype = np.double)
    z = np.zeros(21, dtype = np.double)
    
    for i in range(0,N,1):
        x[i+1]=x[i]-h
    
    #Resolvemos utilizando RK3 de derecha a izquierda
    for i in range(N,0,-1):
        k1 = h*F(x[i],y[i],z[i])
        l1 = h*G(x[i],y[i],z[i])
        
        k2=h*F(x[i]+h/2,y[i]+k1/2,z[i]+l1/2)
        l2=h*G(x[i]+h/2,y[i]+k1/2,z[i]+l1/2)
        
        k3=h*F(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2)
        l3=h*G(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2)
        
        y[i-1]=y[i]+(k1+4*k2+k3)/6
        z[i-1]=z[i]+(l1+4*l2+l3)/6
        
    
    if(imprimir == 1):
        v = np.zeros_like(x)

        for i in range(0,N,1):
            v[i] = np.abs(solucionExacta(x[i])-y[i])

        df = pd.DataFrame({'x_i':x,'u_i':y,'Error':v}, 
                          columns = ['x_i','u_i','Error'])

        df["x_i"]= df["x_i"].apply(lambda x: '{: 0.3f}'.format(x))
        df["u_i"]= df["u_i"].apply(lambda x: '{: 0.7E}'.format(x))
        df["Error"]= df["Error"].apply(lambda x: '{: 0.7E}'.format(x))
        print(df)
        
        fig, axes = plt.subplots(1,1,sharey=True)
        axes.set(xlabel = "$x$", ylabel = "$u(x)$",
                 title = "Método 1, $s = 0.96842632$")
        axes.grid()
        axes.plot(x, y,'o')

        
        return df
    
    else:    
        return y[0]

#Ahora procedemos a aplicar el Método 1

e =  0.000000000000001
i, n = 1, 20

u = np.zeros(200, dtype = np.double)
s = np.zeros(200, dtype = np.double)

s[0], s[1] = 1, 3

u[0] = P1s(s[0],0)

if np.abs(u[0]-1)<e:
    print("s es",s[0],"\n")
    
else:
    while(True):
        u[i] = P1s(s[i],0)
        if np.abs(u[i]-1)<e:
            print("s es",s[i],"\n")
            sf = s[i]
            break
        else:
            s[i+1]=s[i]-(u[i]-1)*((s[i]-s[i-1])/(u[i]-u[i-1]))
        i = i+1

df = P1s(sf,1)

#Enviamos los resultados a un excel llamado ME1E1
df.to_excel("M1E3.xlsx",index = False)

#Tambien se puede generar la tabla automáticamente en TeX
df.to_latex("M1E3.tex",index = False)
