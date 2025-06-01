import ntpath
import scipy.integrate as spi
import scipy as s
import numpy as np
import matplotlib.pyplot as plt

#conditions initiales
k=1
m=1
g=10
z0=2

#resolution de l'edo
# def f(z,t):
#     return((-k/m)*z-g-z0*(k/m))
# z=np.linspace(-5,5,5000)
# Z=spi.odeint(f,z0,z)
# plt.plot(z,Z)
# plt.show()

# #euler explicite
# t=np.linspace(0,5,1000)
# h=t[1]
# z=np.zeros(len(t))
# z[0]=z0
# for i in range(len(t)-1):
#     z[i+1]=z[i]+h*(-(k/m)*z[i])
# plt.plot(t,z)
# plt.show()

# t=np.linspace(0,5,1000)
# h=t[1]
# z=np.zeros(len(t))
# z[0]=z0
# for i in range(len(t)-1):
#     z[i+1]=z[i]+h*(-(k/m)*z[i])
# plt.plot(t,z)
# plt.show()

#conditions initiales
k=1
m=1
g=10
x0=2
v0=0.5
t=np.linspace(0,5,1000)
h=t[1]
X=np.zeros((len(t),2))
X[0,0]=x0
X[0,1]=v0
A=np.array([[0,1],[-1,0]])
for i in range(len(t)-1):
    X[i+1]=X[i]+h*X[i]

plt.plot(t,X)
plt.show()

#autre redaction
# k=1
# m=1   
# g=10
# x0=2
# v0=4
# N=5000 
# X=np.zeros((2,ntpath))
# t=np.linspace(0,5,1000)
# A=np.array([[0,1],[-1,0]])
# for n,tc in enumerate(t[0:-2]):
#     X(0,n+1)
