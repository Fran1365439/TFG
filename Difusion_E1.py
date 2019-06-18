# -*- coding: utf-8 -*-
"""
Created on Sat May 25 11:53:29 2019

@author: Fran
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 12:15:28 2019

@author: Fran
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

def F(x,y,z):
    return z

def G(x,y,z,k, sn_ant, un_ant, xn_ant):
    q = 1/k
    return q*y-g(x, sn_ant, un_ant, xn_ant, q)

def f(x):
    return -1

def g(x, sn_ant, un_ant, xn_ant, q):
    if(np.sum(np.abs(xn_ant)) == 0):
        return f(x)
    
    f2 = interp1d(xn_ant, un_ant, kind='quadratic')

    if x < sn_ant:
        return f2(x)*q+f(x)
    elif x >= sn_ant:
        return f(x)

def P1s(N, s, k, sn_ant, un_ant, xn_ant, imprimir):
    """
    Resolvemos el problema P1s con el método de RK3 de derecha a izquierda.
    Denotamos a la solución por u
    """
    a, b = s, 0
    h = (b-a)/N
    
    x = np.zeros(N+1, dtype = np.double)
    y = np.zeros(N+1, dtype = np.double)
    z = np.zeros(N+1, dtype = np.double)
    
    for i in range(0,N,1):
        x[i+1]=x[i]-h
    
    #Resolvemos utilizando RK3 de derecha a izquierda
    for i in range(N,0,-1):
        k1 = h*F(x[i],y[i],z[i])
        l1 = h*G(x[i],y[i],z[i],k, sn_ant, un_ant, xn_ant)
        
        k2=h*F(x[i]+h/2,y[i]+k1/2,z[i]+l1/2)
        l2=h*G(x[i]+h/2,y[i]+k1/2,z[i]+l1/2,k, sn_ant, un_ant, xn_ant)
        
        k3=h*F(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2)
        l3=h*G(x[i]+h,y[i]-k1+2*k2,z[i]-l1+2*l2,k, sn_ant, un_ant, xn_ant)
        
        y[i-1]=y[i]+(k1+4*k2+k3)/6
        z[i-1]=z[i]+(l1+4*l2+l3)/6
          
    if(imprimir == 2):
        return y
    if(imprimir == 3):
        return x        
    else:    
        return y[0]

#Ahora procedemos a aplicar el Método 1
def MethodI(N, s1, s2, U, k, sn_ant, un_ant, xn_ant, ret):
    e =  0.00000001
    i = 1
    
    u = np.zeros(200, dtype = np.double)
    s = np.zeros(200, dtype = np.double)
    
    s[0], s[1] = s1, s2
    
    u[0] = P1s(N, s[0], k, sn_ant, un_ant, xn_ant, 0)
    
    if np.abs(u[0]-U)<e:
#        print("s es",s[0],"\n")
        sf = s[0]
        
    else:
        while(True):
            u[i] = P1s(N, s[i], k, sn_ant, un_ant, xn_ant, 0)
            if np.abs(u[i]-U)<e:
#                print("s es",s[i],"\n")
                sf = s[i]
                break
            else:
                s[i+1]=s[i]-(u[i]-U)*((s[i]-s[i-1])/(u[i]-u[i-1]))
            i = i+1
    
#    df = P1s(sf,1)
    
    if(ret == 0):
        return P1s(N, sf, k, sn_ant, un_ant, xn_ant, 2) #Return valores un
    elif(ret == 1):
        return sf #Return valor de sn
    elif(ret == 2):
        return P1s(N, sf, k, sn_ant, un_ant, xn_ant, 3) #Return valores xn
    
    
    
def ODEuij(N, j, b, xn, un, sn, k):
    s1, s2 = b, b + 0.3
    
    xn[j+1,:] = MethodI(N = N, s1 = s1, s2 = s2, U = 1/2, k = k, sn_ant = sn[j], un_ant = un[j,:], xn_ant = xn[j,:], ret = 2)
    un[j+1,:] = MethodI(N = N, s1 = s1, s2 = s2, U = 1/2, k = k, sn_ant = sn[j], un_ant = un[j,:], xn_ant = xn[j,:], ret = 0)
    sn[j+1] =   MethodI(N = N, s1 = s1, s2 = s2, U = 1/2, k = k, sn_ant = sn[j], un_ant = un[j,:], xn_ant = xn[j,:], ret = 1)
    
    b = sn[j+1]
    
    for i in range(0,N,1):
        if(xn[j+1,i] < sn[j+1]):
            un[j+1,i] = un[j+1,i]
            
        elif(xn[j+1,i] >= sn[j+1]):
            un[j+1,i] = 0


#    df = pd.DataFrame({'x_i':xn[j+1,:],'u_i':un[j+1,:]}, 
#                          columns = ['x_i','u_i'])
#
#    df["x_i"]= df["x_i"].apply(lambda x: '{: 0.3f}'.format(x))
#    df["u_i"]= df["u_i"].apply(lambda x: '{: 0.7E}'.format(x))
#    print(df)
    
#    fig, axes = plt.subplots(1,1,sharey=True)
#    axes.set(xlabel = "$x$", ylabel = "$u(x)$",
#             title = " ")
#    axes.grid()
#    axes.plot(xn[j+1,:], un[j+1,:],'o') 
    
#    print("xn \n", xn)
#    print("un \n", un)
#    print("sn \n", sn)



N = 11
N1 = 40
k = 0.02
T = 0.42

xn = np.zeros((N1,N+1))
un = np.zeros((N1,N+1)) #matriz que contiene las soluciones para cada n. i = iteracion, j = valor en cada punto
sn =  np.zeros(N1)

un[0,:] = 0 #primer valor es 0
b = 1

for j in range(0,N1-1,1):
    ODEuij(N, j, b, xn, un, sn, k)
    
    if k*(j+1) == T:
        break
  
    
#plots soluciones en
        
#t = 0.02
#######
df1 = pd.DataFrame({'x_i':xn[1,:],'u_i':un[1,:]}, 
                      columns = ['x_i','u_i'])

df1["x_i"]= df1["x_i"].apply(lambda x: '{: 0.3f}'.format(x))
df1["u_i"]= df1["u_i"].apply(lambda x: '{: 0.7E}'.format(x))

fig, axes = plt.subplots(1,1,sharey=True)
axes.set(xlabel = "$x$", ylabel = "$u(x)$",
         title = " ")
axes.grid()
axes.plot(xn[1,:], un[1,:],'o')
print(df1)


#t = 0.1
########
df2 = pd.DataFrame({'x_i':xn[6,:],'u_i':un[6,:]}, 
                      columns = ['x_i','u_i'])

df2["x_i"]= df2["x_i"].apply(lambda x: '{: 0.3f}'.format(x))
df2["u_i"]= df2["u_i"].apply(lambda x: '{: 0.7E}'.format(x))
print(df2)

fig, axes = plt.subplots(1,1,sharey=True)
axes.set(xlabel = "$x$", ylabel = "$u(x)$",
         title = " ")
axes.grid()
axes.plot(xn[6,:], un[6,:],'o')


#t = 0.4
########
df3 = pd.DataFrame({'x_i':xn[21,:],'u_i':un[21,:]}, 
                      columns = ['x_i','u_i'])

df3["x_i"]= df3["x_i"].apply(lambda x: '{: 0.3f}'.format(x))
df3["u_i"]= df3["u_i"].apply(lambda x: '{: 0.7E}'.format(x))
print(df3)

fig, axes = plt.subplots(1,1,sharey=True)
axes.set(xlabel = "$x$", ylabel = "$u(x)$",
         title = " ")
axes.grid()
axes.plot(xn[21,:], un[21,:],'o')

#plot s(t)
time = np.zeros(22)

for i in range(0,22,1):
    time[i] = 0.02*i

time

fig, axes = plt.subplots(1,1,sharey=True)
axes.set(xlabel = "$t$", ylabel = "$s(t)$",
         title = " ")
axes.grid()
axes.plot(time, sn[0:22],'o')


xss = np.array([0.000,
               0.091,
               0.182,
               0.273,
               0.364,
               0.455,
               0.545,
               0.636,
               0.727,
               0.818,
               0.909,
               1.000])
uss = np.array([5.0000000E-01,
                4.1322314E-01,
                3.3471074E-01,
                2.6446281E-01,
                2.0247934E-01,
                1.4876033E-01,
                1.0330579E-01,
                6.6115702E-02,
                3.7190083E-02,
                1.6528926E-02,
                4.1322314E-03,
                0.0000000E+00]) 

fig, axes = plt.subplots(2,2,sharey=True, sharex=True)
#plt.text(-12, 12.5, "Euler Explícito C.F. Dirichlet",
#         horizontalalignment='center',
#         fontsize=18)
plt.subplots_adjust(hspace = 0.55)
axes[0,0].set_title("T = 0.02")
axes[0,0].set(xlabel = "$x$", ylabel = "$u(x)$")
axes[0,0].plot(xn[1,:], un[1,:],'o')
axes[0,1].set_title("T = 0.1")
axes[0,1].set(xlabel = "$x$", ylabel = "$u(x)$")
axes[0,1].plot(xn[6,:], un[6,:],'o')
axes[1,0].set_title("T = 0.42")
axes[1,0].set(xlabel = "$x$", ylabel = "$u(x)$")
axes[1,0].plot(xn[21,:], un[21,:],'o')
axes[1,1].set_title("Steady-State")
axes[1,1].set(xlabel = "$x$", ylabel = "$u(x)$")
axes[1,1].plot(xss, uss,'o')

#Tambien se puede generar la tabla automáticamente en TeX
#df.to_latex("DIF_T01_E1.tex",index = False)


