import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

M = np.array((2,2))
I = np.eye(2)

# search of the critical point
# with consider A = 1.0  and make variation of B
A = 1.0
n_plt = 40+1

B_plt = np.linspace(0,3,n_plt) #  we choose 0< B <3   

lbd_plt = np.zeros(n_plt) # i'm looking for the max real part eigenvalues

def GetEigen(A,B):
    M=[[B-1,A**2],[-B,-A**2]] # the jacobian at the stability point
    lbd = sp.symbols('lbd') # to use symbolic computation lbd is a symbal 
    DG=sp.Matrix(M)-lbd*sp.Matrix(I) # M - lbd I 
    eqn_char = DG.det() # its dterminant
    eigs=sp.solve(eqn_char,lbd) # its eigenvalues numerical from symbolic
    g = sp.lambdify((), eigs, modules='numpy') # convert to numpy
    max_eig=max(np.real(g())) # get the maximum real part 
    return max_eig

nA,nB=20,20
A = np.linspace(0,10,nA)
B = np.linspace(0,10,nB)
Av, Bv = np.meshgrid(A, B, indexing='ij')
stability=np.chararray(np.shape(Av))
x,y,m=[],[],[]
for a in A:
    for b in B:
        x.append(a)
        y.append(b)
        res=GetEigen(a,b)
        if (res>0) :
            m='+'
        else :
            m='x'
        
        plt.scatter(a,b, marker=m,c=0)

plt.show()









