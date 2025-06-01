import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as so
import sympy as sym

#Partie 2 question 5

B0=4
A0=10
x0=A0
y0=B0/A0

n=100
t=np.linspace(0,10,n)
h=t[1]-t[0]
V=np.zeros((n,2))
V[0,0]=x0
V[0,1]=y0

for i in range(n-1):
    V[i+1,0]= ((-1-B0)*V[i,0]+((V[i,0])**2)*V[i,1])*h + A0
    V[i+1,1]= (B0*V[i,0]-((V[i,0])**2)*V[i,1])*h

for i in range (2):
    plt.figure(i)
    plt.plot(t,V[:,i])
    plt.show


#Partie 3 question 2

# x0=0
# y0=0
# z0=0
# v=[0.9,1.3,1.52]
# n=100
# t=np.linspace(0,10,n)
# h=t[1]-t[0]
# V=np.zeros((n,3))
# V[0,0]=x0
# V[0,1]=y0
# V[0,2]=z0
# for j in range(3):
#     for i in range(n):
#         V[i+1,0]= 
#         V[i+1,1]=  
#         V[i+1,2]=       
