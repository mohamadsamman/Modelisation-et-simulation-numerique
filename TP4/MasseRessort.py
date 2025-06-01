import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

T=2*np.pi*10# time of the simulation
N=1600*2 # number of subdivisation
h=T/N
t=np.linspace(0,T,N) # definition of the array containing time-step
q,p=np.zeros(N+1),np.zeros(N+1) #

# Euler explicite
q[0],p[0]=0,1 # (position,vitesse)
for n,tc in enumerate(t):
    q[n+1] = q[n] + h*p[n]
    p[n+1] = p[n] - h*q[n]

f = open("MasseRessort_EE.dat","w")
em0= 0.5*p[0]*p[0] + 0.5*q[0]*q[0]
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n] - em0
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em))
f.close()

    
# Euler implicite
q[0],p[0]=0,1
for n,tc in enumerate(t):
    q[n+1] = 1/(1+h**2)*(   q[n] + h*p[n])
    p[n+1] = 1/(1+h**2)*(-h*q[n] +   p[n])


f = open("MasseRessort_EI.dat","w")
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n]
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em))
f.close()

q[0],p[0]=0,1
for n,tc in enumerate(t):
    Minv= inv( np.array([[1.,-h/2],[h/2,1.]] ) )
    dg= np.array([q[n] + 0.5*h*p[n],p[n] - 0.5*h*q[n]])
    res= np.matmul(Minv,dg)
    q[n+1] = res[0]
    p[n+1] = res[1]

f = open("MasseRessort_CN.dat","w")
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n]
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em))
f.close()

q[0],p[0]=0,1
for n,tc in enumerate(t):
    Minv= inv( np.array([[1.,0.],[h,1.]] ) )
    dg= np.array([q[n] + h*p[n],p[n] ])
    res= np.matmul(Minv,dg)
    q[n+1] = res[0]
    p[n+1] = res[1]

f = open("MasseRessort_Symplectique.dat","w")
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n] + h*p[n]*q[n]*0.5
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em))
f.close()



q[0],p[0]=0,1
for n,tc in enumerate(t):
    q[n+1] = q[n] + h*p[n]
    p[n+1] = p[n] - h*q[n+1]

f = open("MasseRessort_ESa.dat","w")
em0= 0.5*p[0]*p[0] + 0.5*q[0]*q[0]
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n] #+ h*p[n]*q[n]*0.5
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em-em0))
f.close()


q[0],p[0]=0,1
for n,tc in enumerate(t):
    p[n+1] = p[n] - h*q[n]
    q[n+1] = q[n] + h*p[n+1]


f = open("MasseRessort_ESb.dat","w")
em0= 0.5*p[0]*p[0] + 0.5*q[0]*q[0]
for n,tc in enumerate(t):
    em = 0.5*p[n]*p[n] + 0.5*q[n]*q[n] #+ h*p[n]*q[n]*0.5
    f.write("%f %e %e %e\n"%(tc,q[n],p[n],em-em0))
f.close()

em =  0.5*p*p + 0.5*q*q - em0
print(np.max(em))

exit()




