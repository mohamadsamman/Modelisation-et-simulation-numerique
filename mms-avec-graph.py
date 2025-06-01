#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np
#import tecplot_writer as tp

# domaine de calcul
Lx,Ly=2*np.pi,2*np.pi
Nx,Ny=40,40
hx,hy=Lx/Nx,Ly/Ny

# maillage face et centre en x
Xf=np.linspace(0,Lx,Nx, endpoint=False)
Xc=np.linspace(hx/2.,Lx-hx/2.,Nx)

# maillage face et centre en y
Yf=np.linspace(0,Ly,Ny, endpoint=False)
Yc=np.linspace(hy/2.,Ly-hy/2.,Ny)

# declaration des variables comme champs de données indicés comme (i,j)
size=(Nx,Ny)
u,v,pres,fi=np.zeros(size),np.zeros(size),np.zeros(size),np.zeros(size)
u_star,v_star=np.zeros(size),np.zeros(size)
NLu,NLv,Lu,Lv=np.zeros(size),np.zeros(size),np.zeros(size),np.zeros(size)
sfi=np.zeros(size)

# declaration des variables comme vecteur (p=i+j*nx) pour résoudre la correction de pression
b,x=np.zeros(Nx*Ny),np.zeros(Nx*Ny) 

# pour traiter les conditions aux limites périodiques
II = np.arange(Nx)
IP = np.arange(Nx)+1
IM = np.arange(Nx)-1
IP[Nx-1],IM[0]=0,Nx-1
#print(IP)
#print(IM)
JJ = np.arange(Ny)
JP = np.arange(Ny)+1
JM = np.arange(Ny)-1
JP[Ny-1],JM[0]=0,Ny-1


# calcul des termes non-linéaires NLu & NLv
def get_NonLinear_terms(u,v):
    NLu , NLv = np.zeros(np.shape(u)) , np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im,ip,jm,jp=IM[i],IP[i],JM[j],JP[j]
            ul = 0.5*( u[im,j] + u[i,j] )
            ur = 0.5*( u[ip,j] + u[i,j] )
          
            ub = 0.5*( u[i,jm] + u[i,j] )
            ut = 0.5*( u[i,jp] + u[i,j] )
          
            up = u[i,j]
            vp = 0.25*( v[i,j] + v[i,jp] + v[im,jp] + v[im,j] )
          
            NLu[i,j] = up*(ur-ul)/hx  + vp*(ut-ub)/ hy

            vb = 0.5*( v[i,jm] + v[i,j] )
            vt = 0.5*( v[i,jp] + v[i,j] )
            
            vl = 0.5*( v[im,j] + v[i,j] )
            vr = 0.5*( v[ip,j] + v[i,j] )
            
            vp = v[i,j]
            up = 0.25*( u[i,j] + u[ip,j] + u[i,jm] + u[ip,jm] )
            
            NLv[i,j] = up*(vr-vl)/hx  + vp*(vt-vb)/hy
    
    return [NLu,NLv]

# calcul des termes non-linéaires Lu & Lv
def get_Linear_terms(u,v):
    Lu , Lv = np.zeros(np.shape(u)) , np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im,ip,jm,jp=IM[i],IP[i],JM[j],JP[j]
            Lu[i,j] =           ( u[ip,j] - 2*u[i,j] + u[im,j] )/hx**2
            Lu[i,j] = Lu[i,j] + ( u[i,jp] - 2*u[i,j] + u[i,jm] )/hy**2
            
            Lv[i,j] =           ( v[ip,j] - 2*v[i,j] + v[im,j] )/hx**2
            Lv[i,j] = Lv[i,j] + ( v[i,jp] - 2*v[i,j] + v[i,jm] )/hy**2

    return [Lu,Lv]

# calcul de la divergence
def get_divergence(u,v):
    Dx , Dy = np.zeros(np.shape(u)) , np.zeros(np.shape(v))
    for i in II:
        for j in JJ:
            im,ip,jm,jp=IM[i],IP[i],JM[j],JP[j]
            Dx[i,j] = ( u[ip,j] - u[i,j] )/hx
            Dy[i,j] = ( v[i,jp] - v[i,j] )/hy
    div = Dx+Dy
    return div

def get_gradient(fi):

    Gx , Gy = np.zeros(np.shape(fi)) , np.zeros(np.shape(fi))
    for i in II:
        for j in JJ:
            im,ip,jm,jp=IM[i],IP[i],JM[j],JP[j]
            Gx[i,j] = ( fi[i,j] - fi[im,j] )/hx
            Gy[i,j] = ( fi[i,j] - fi[i,jm] )/hy
    return [Gx,Gy]


def index(i,j):
    p=i+j*Nx
    return p
    
def BuildPoisson(hx,hy,nx,ny):
    M=np.zeros((nx*ny,nx*ny))
    for i in range(nx):
        for j in range(ny):
            im,ip,ii=IM[i],IP[i],II[i]
            jm,jp,jj=JM[j],JP[j],JJ[j]
            row = index(ii,jj)
            M[row, index(ii,jj)] = -2/hx**2-2/hy**2
            M[row, index(im,jj)] = 1/hx**2
            M[row, index(ip,jj)] = 1/hx**2
            M[row, index(ii,jm)] = 1/hy**2
            M[row, index(ii,jp)] = 1/hy**2

    M[0,:]=M[0,:]*0
    M[0,0]=1
    
    return M

# en précalcul j'assemble le système linéaire M
M=BuildPoisson(hx,hy,Nx,Ny)
# en précalcul on calcule l'inverse de M
Minv=np.linalg.inv(M)
#print(np.matmul(Minv, M))


nu = 1e-1 # viscosité
dt = 0.25*hx**2/nu # stabilité
#
xc,yc= np.meshgrid(Xc,Yc, indexing='ij')
xu,yu= np.meshgrid(Xf,Yc, indexing='ij')
xv,yv= np.meshgrid(Xc,Yf, indexing='ij')

#from Ns_inc_2d import *

from numpy import cos as Cos
from numpy import sin as Sin

a=2*np.pi #problème instationnaire
#a = 0. 
def set_su(t,x,y):
    su=Cos(x)*( -Cos(a*t)**2*Sin(x) - a*Sin(a*t)*Sin(y) + 2*Cos(a*t)*(2*Sin(x) + nu*Sin(y)))
    return su

def set_sv(t,x,y):
  sv=Cos(y)*((-2*nu*Cos(a*t) + a*Sin(a*t))*Sin(x) - (-4 + Cos(a*t))*Cos(a*t)*Sin(y))
  return sv

def set_u(t,x,y):
  u=+Cos(x)*Sin(y)*Cos(a*t)
  return u

def set_v(t,x,y):
  v=-Sin(x)*Cos(y)*Cos(a*t)
  return v

def set_p(t,x,y):
  pres=-(Cos(2*x)+Cos(2*y))*Cos(a*t)
  return pres

  



t,Tmax=0.,40.
# initialisation MMS
u = set_u(t,xu,yu)
v = set_v(t,xv,yv)
pres = set_p(t,xc,yc)

f = open("erreur.dat","w")
while(t<Tmax):
    t=t+dt
    NLu,NLv=get_NonLinear_terms(u,v)

    Lu,Lv = get_Linear_terms(u,v)

    dpdx,dpdy=get_gradient(pres)

    Su=set_su(t,xu,yu)
    Sv=set_sv(t,xv,yv)

    u_star = u + dt*( - dpdx - NLu + nu*Lu + Su )
    v_star = v + dt*( - dpdy - NLv + nu*Lv + Sv )

    div=get_divergence(u_star,v_star)

    div=+div/dt

    for i in II: 
        for j in JJ:
            row = index(i,j)
            b[row] = div[i,j] # passage de la notation tableau (i,j) en ve vecteur p = i+j*nx  

    b[0]=0
    # très couteux :-(
    x = Minv @ b  # resolution  de type M^{-1} b      

    for i in II:
        for j in JJ:
            row = index(i,j)  # passage de la notation  vecteur p = i+j*nx en tableau (i,j) 
            fi[i,j] = x[row]

    Gx , Gy = get_gradient(fi)
    u = u_star - dt*Gx
    v = v_star - dt*Gy
    pres = pres + fi

    div = get_divergence(u,v)


    ue=set_u(t,xu,yu)
    ve=set_v(t,xv,yv)
    pe=set_p(t,xc,yc)

    err_u = np.max(np.abs(ue-u))
    err_v = np.max(np.abs(ve-v))
    div_max =  np.max(np.abs(div))
    #print("%f %f %f %f %e"%(t,hx,err_u,err_v,div_max))
    #f.write("%f %f %f %f %e\n"%(t,hx,err_u,err_v,div_max))



# pour tracer une iso
import matplotlib.pyplot as plt

plt.figure(1)
plt.contourf(xu,yu,u)
plt.title('iso u')
plt.show()
plt.figure(2)
plt.contourf(xv,yv,v)
plt.title('iso v')
plt.show()
# pour plotter avec tecplot
#tp.tecplot_writer('test2d.dat', {'u': u,'v':v,'pres':fi}, Xf, Yf)




omega=np.copy(u)
for i in II: 
    for j in JJ:
        omega[i,j]= (v[i,j]-v[IM[i],j])/hx-(u[i,j]-u[i,JM[j]])/hy
plt.figure(3)
plt.contourf(Xf,Yf,omega)
plt.title('Vorticité')
plt.show()

ue,ve=np.cos(xu)*np.sin(yu),-np.sin(xv)*np.cos(yv)
print(np.max(np.abs(ue-u)))
err=ue-u

#matplotlib.pyplot.contourf(*args, data=None, **kwargs)


# In[ ]:





# In[ ]:




