import sympy as sp
import numpy as np
import matplotlib.pyplot as mp
import time 

# we consider

def euler_explicite(F,ti,tf,x0,N):
    h=(tf-ti)/N                # time step
    t,x=np.linspace(ti,tf,N+1),np.zeros((2,N+1)) # allocate time subdivision and solution
    x[0:1,0]=x0[0:1] # initial condition
    for n,tc in enumerate(t[:-1]):
        f=F(t,x[:,n])
        x[0,n+1] = x[0,n] + h*f[0]
        x[1,n+1] = x[1,n] + h*f[1]
    return [t,x]

def euler_implicite(F,DF,ti,tf,x0,N,tol=1e-6,nl_it_max=50):
    h=(tf-ti)/N                 # pas d'espace
    t,x=np.linspace(ti,tf,N+1),np.zeros((2,N+1))  # le temps et la solution vide
    x[0:1,0]=x0[0:1] # initial condition
    Id,rhs = np.eye(2),np.array([0.,0.])
    
    for n,tc in enumerate(t[:-1]):
        xnew,xn = x[:,n],x[:,n] #
        nl_it = 0
        rhs = (xnew-xn)-h*F(t[n]+h,xnew)
        while ((nl_it<nl_it_max) and  np.max(np.abs(rhs))>tol):
            rhs = (xnew-xn)-h*F(t[n]+h,xnew) 
            jac = Id - h*DF(t[n]+h,xnew)
            inc = -np.linalg.solve(jac,rhs)
            xnew = xnew + inc
            nl_it = nl_it + 1
        x[:,n+1]=xnew

    return t,x





A , B = 1.5 , 1.5**2 + 1 - 1 


def F(t,X) :
    x , y = X[0] , X[1]
    Fx = A - (B+1)*x + x**2*y
    Fy =     (B  )*x - x**2*y
    res = np.array([Fx,Fy])
    return res

def DF(t,X) :
    x , y = X[0] , X[1]
    Fxx,Fxy = -(B+1) + 2*x*y , +x**2
    Fyx,Fyy =   B-2*x*y      , -x**2
    res = np.array([[Fxx,Fxy],[Fyx,Fyy]])
    return res

x0=np.array([A,B/A])
x0=x0+1e-4


#M=np.array([[B-1,A**2],[-B,-A**2]])
#lbd,V = np.linalg.eig(M)
#print(lbd)



tic = time.perf_counter()
t,x=euler_explicite(F,0.,100.,x0,1000)
toc = time.perf_counter()
print(f"euler_explicite {toc - tic:0.4f} seconds")

f = open("euler_explicite_1.dat","w")
for n,tc in enumerate(t):
    f.write("%f %e %e \n"%(tc,x[0,n],x[1,n]))
f.close()

tic = time.perf_counter()
t,x=euler_implicite(F,DF,0.,100.,x0,1000)
toc = time.perf_counter()
print(f"euler_implicite {toc - tic:0.4f} seconds")

f = open("euler_implicite_1.dat","w")
for n,tc in enumerate(t):
    f.write("%f %e %e \n"%(tc,x[0,n],x[1,n]))
f.close()


