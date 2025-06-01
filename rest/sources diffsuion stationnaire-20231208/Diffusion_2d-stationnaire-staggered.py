import numpy as np
import matplotlib.pyplot as plt

xmin , xmax , nx = 0. , 1. , 100
ymin , ymax , ny = 0. , 0.1 , 10

hx , hy = (xmax-xmin)/nx , (ymax-ymin)/ny
xc = np.linspace(xmin+.5*hx,xmax-0.5*hx,nx)
yc = np.linspace(ymin+.5*hy,ymax-0.5*hy,ny)

ii,jj=np.arange(nx),np.arange(ny)
nxy=nx*ny
M=np.zeros((nxy,nxy))
b=np.zeros((nxy))
x=np.zeros((nxy))
for j in jj:
    for i in ii:
        row = i + j*(nx)

        b[row] = -1
        
        Ap =  i + j*(nx)

        Aw,Ae = i-1 + j*nx , i+1 + j*nx
        M[row,Ap] = 0
        M[row,Ap] = M[row,Ap] - 2/hx**2
        if (i==0) :
            M[row,Ap] =  M[row,Ap] - 1/hx**2
            M[row,Ae] =  M[row,Ae] + 1/hx**2
        elif (i==(nx-1)):
            M[row,Ap] =  M[row,Ap] - 1/hx**2
            M[row,Aw] =  M[row,Aw] + 1/hx**2
        else:
            M[row,Aw] =  M[row,Aw] + 1/hx**2
            M[row,Ae] =  M[row,Ae] + 1/hx**2

        As,An = i + (j-1)*nx , i + (j+1)*nx

        M[row,Ap] = M[row,Ap] - 2/hy**2
        if (j==0) :
            M[row,Ap] =  M[row,Ap] - 1/hy**2
            M[row,An] =  M[row,An] + 1/hy**2
        elif (j==(ny-1)):
            M[row,Ap] =  M[row,Ap] - 1/hy**2
            M[row,As] =  M[row,As] + 1/hy**2
        else:
            M[row,As] =  M[row,As] + 1/hy**2
            M[row,An] =  M[row,An] + 1/hy**2
print(M)
x=np.linalg.solve(M,b)
fi=np.zeros((nx,ny))
for j in jj:
    for i in ii:
        row = i + j*(nx)
        fi[i,j] = x[row]

xc_v, yc_v = np.meshgrid(xc, yc,indexing='ij')
plt.contourf(xc_v,yc_v,fi)
plt.show()
            
print(fi)
