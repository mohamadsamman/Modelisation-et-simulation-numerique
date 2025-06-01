import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import matplotlib.pyplot as plt
T=10# time of the simulation
N=1600 # number of subdivisation
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
#def A
A=np.zeros((2*n,2*n))
A[:n,n:]=np.eye(n)
A[n:,:n]=w02*V

X = np.zeros((2*n,len(t)))
X[n:,0]= 1
for n in range(len(t)-1):
    X[:,n+1]=X[:,n]+h*np.dot(A,X[:,n])


plt.plot(t,X[n,:])
plt.show()





