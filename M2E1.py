# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 16:40:23 2019

@author: Fran
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def F(x,y,z):
    return z

def G(x,y,z):
    return np.sqrt(x)

def solucionExacta(x):
    """
    Solución exácta del problema evaluada en el punto x.
    El valor de s se ha de definir dentro.
    """
    s = (5/2)**(2/5)
    return (4/15)*(x**2.5-s**2.5)-(2/3)*s**1.5*(x-s)

def P1s(gamma, s, ret):
    """
    Resolvemos el problema P1s con el método de RK3 de izquierda a derecha.
    Denotamos a la solución por u
    Hay tres posibles outputs segun el valor de ret = 0,1,2
    """
    a, b = 0, s
    N = 20
    h = (b-a)/N
    
    x = np.zeros(21, dtype = np.double)
    y = np.zeros(21, dtype = np.double)
    z = np.zeros(21, dtype = np.double)
    
    for i in range(0,N,1):
        x[i+1]=x[i]+h
    
    #Condiciones iniciales
    y[0]=1
    z[0]=gamma

    #Resolvemos utilizando RK3 de derecha a izquierda
    for i in range(0,N,1):
        k1 = h*F(x[i],y[i],z[i])
        l1 = h*G(x[i],y[i],z[i])
        
        k2=h*F(x[i]+h/2,y[i]+k1/2,z[i]+l1/2)
        l2=h*G(x[i]+h/2,y[i]+k1/2,z[i]+l1/2)
        
        k3=h*F(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2)
        l3=h*G(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2)
        
        y[i+1]=y[i]+(k1+4*k2+k3)/6
        z[i+1]=z[i]+(l1+4*l2+l3)/6
        
    
    if(ret == 2):
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
                 title = "Método 2, $s = 1.44269790$")
        axes.grid()
        axes.plot(x, y,'o')

        return df
    
    elif(ret == 0):
#        print(y[N])
        return y[N]
    
    elif(ret == 1):
#        print(z[N])
        return z[N]

#Ahora procedemos a aplicar el Método 1

e =  0.000000000000001
i = 1

u = np.zeros(200, dtype = np.double)
u1 = np.zeros(200, dtype = np.double)
s = np.zeros(200, dtype = np.double)
d = np.zeros(200, dtype = np.double)

Gnn = np.zeros(2)
Gn1n = np.zeros(2)
Gnn1 = np.zeros(2)
J = np.zeros((2,2))
Hd = np.zeros((2,2))
Hs = np.zeros((2,2))

s[0], s[1] = 1, 3
d[0], d[1]= -1.5, -3.7

u[0] = P1s(d[0],s[0],0)
u1[0]= P1s(d[0],s[0],1)

if ((np.abs(u[0])+np.abs(u1[0]))<e):
    print("s es",s[0],"\n")
    
else:
    while(True):     
        u[i] = P1s(d[i],s[i],0)
        u1[i]= P1s(d[i],s[i],1)

        
        if ((np.abs(u[i])+np.abs(u1[i]))<e):
            print("s es",s[i],"\n")
            print("d es",d[i],"\n")
            sf = s[i]
            df = d[i]
            break
            
        else:
            u[i] = P1s(d[i-1],s[i],0)
            u1[i]= P1s(d[i-1],s[i],1)
            
            if ((np.abs(u[i])+np.abs(u1[i]))<e):
                print("s es",s[i],"\n")
                print("d es",d[i],"\n")
                sf = s[i]
                df = d[i]
                break
                
            else:
                u[i] = P1s(d[i],s[i-1],0)
                u1[i]= P1s(d[i],s[i-1],1)
                
                if ((np.abs(u[i])+np.abs(u1[i]))<e):
                    print("s es",s[i],"\n")
                    print("d es",d[i],"\n")
                    sf = s[i]
                    df = d[i]
                    break
                
                else:
                    Gnn[0] = P1s(d[i],s[i],0)
                    Gnn[1]= P1s(d[i],s[i],1)
                    
                    Gn1n[0] = P1s(d[i-1],s[i],0)
                    Gn1n[1]= P1s(d[i-1],s[i],1)
                    
                    Gnn1[0] = P1s(d[i],s[i-1],0)
                    Gnn1[1]= P1s(d[i],s[i-1],1)
                    
                    J[0][0]=(Gnn[0]-Gn1n[0])/(d[i]-d[i-1]);
                    J[0][1]=(Gnn[0]-Gnn1[0])/(s[i]-s[i-1]);
                    J[1][0]=(Gnn[1]-Gn1n[1])/(d[i]-d[i-1]);
                    J[1][1]=(Gnn[1]-Gnn1[1])/(s[i]-s[i-1]);
                    
                    
                    Hd[0][0]=Gnn[0];
                    Hd[0][1]=(Gnn[0]-Gnn1[0])/(s[i]-s[i-1]);
                    Hd[1][0]=Gnn[1];
                    Hd[1][1]=(Gnn[1]-Gnn1[1])/(s[i]-s[i-1]);
                    Hs[0][0]=(Gnn[0]-Gn1n[0])/(d[i]-d[i-1]);
                    Hs[0][1]=Gnn[0];
                    Hs[1][0]=(Gnn[1]-Gn1n[1])/(d[i]-d[i-1]);
                    Hs[1][1]=Gnn[1];
                    
                    
                    d[i+1]=d[i]-(Hd[0][0]*Hd[1][1]-Hd[1][0]*Hd[0][1])/(J[0][0]*J[1][1]-J[1][0]*J[0][1]);
                    s[i+1]=s[i]-(Hs[0][0]*Hs[1][1]-Hs[1][0]*Hs[0][1])/(J[0][0]*J[1][1]-J[1][0]*J[0][1]);
                    i = i+1

dfT = P1s(df,sf,2)

#Enviamos los resultados a un excel llamado ME1E1
dfT.to_excel("M2E1.xlsx",index = False)

#Tambien se puede generar la tabla automáticamente en TeX
dfT.to_latex("M2E1.tex",index = False)