import numpy as np  
import matplotlib.pyplot as p
import math as m


#euler explicite
# def f(t):
#     return((1/2)*(t**2))

# for i in range(len(x)-1):
#     y[i+1]=y[i]+x[i]*y[i]

# for i in range(len(x)-1):
#     y[i+1]=f(x[i])

# vec_emax=np.zeros(100)
# for k in range(100):
#     x=np.linspace(0,10+k,k)
#     x0=0
#     y=1*x
#     y[0]=0
#     yex=1*1*x

#     for i in range(len(x)-1):
#         y[i+1]=y[i]+x[i]*y[i]

#     for i in range(len(x)-1):  
#         y[i+1]=f(x[i+1])
#         emax=0
#         for i in range(len(y)):
#             e=abs(y[i]-yex[i])
#             if (e>emax):
#                 emax=e
#A FINIR


#heunn
def t3(t):
    return(t**3)
def f(t):
    return(1-m.exp(-t))
def F(x):
    return(1-x)

 
h=np.arange(0.01,0.51,0.01)
k=len(h)
E=np.ones((1,k))
H=1*E
print(E)
for j in range(len(h)):
    x= np.arange(0,10+h[j],h[j])
    y= 1*x
    yex=14*x
    y[0]=0
    yex[0]=0
    for i in range(len(x)-1):
        yex[i+1]=f(x[i+1])

    for i in range(len(x)-1):
        y[i+1]=y[i]+((h[j])/2)*(F(y[i])+y[i]+h[j]*F(y[i])) #ERREUR ICI
    
    for i in range(len(x)-1):  
        emax=0
        for i in range(len(y)):
            e=abs(y[i]-yex[i])
            if (e>emax):
                emax=e
    E[0,j]=emax
    H[0,j]=h[j]
for i in range(len(E)):
    E[0,i]=np.log(E[0,i])
    H[0,i]=np.log(H[0,i])
p.plot(H,E)
p.show()