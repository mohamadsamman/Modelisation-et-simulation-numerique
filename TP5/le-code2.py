import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import matplotlib.pyplot as plt
T=10# time of the simulation
N=1000 # number of subdivisation
h=T/N
t=np.linspace(0,T,N)
#dim n
n=4
#def de V
V=2*np.eye(n)-np.diag(np.ones(n-1),k=1)-np.diag(np.ones(n-1),k=-1)
V[0,n-1]=1
V[n-1,0]=1
#w0²
w02=1
c=0
tc=0.01
p,q=np.zeros(n),np.zeros(n)
while (tc<T and c<1000):
    c=c+1
    tc=tc+0.1
    for i in range(n):
        q[i+1]=q[i]+0.1*p[i]
    for i in range (2,n):
        p[i+1]=p[i]-0.1*(-q[i-1]+2*q[i]-q[i+1])
    p[0]-0.1*(2*q[1]-q[2])
    p[n-1]-0.1*(2*q[n-2]-q[n-1])