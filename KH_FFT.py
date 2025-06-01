import numpy as np
import tecplot_writer as tp

# domaine de calcul
Lx,Ly=2.,1.
Nx,Ny=60,60
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

JJ = np.arange(Ny)
JP = np.arange(Ny)+1
JM = np.arange(Ny)-1
JP[Ny-1],JM[0]=0,Ny-1

def init_field():
  u , v = np.zeros(size) , np.zeros(size)
  Pj,Rj,Ax = 20 , 0.25 , 0.5
  lbd_x = 0.5*Lx
  for i in II:
      for j in JJ:
        u[i,j] = 0.5*( 1 + np.tanh( 0.5*Pj*( 1 - np.abs(Ly/2. - Yc[j])/Rj ) ) )
        u[i,j] = u[i,j]*( 1 + Ax*np.sin(2*np.pi*Xf[i]/lbd_x) )
  return[u,v]

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

# calcul des termes linéaires Lu & Lv
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


def BuildTridiag(hx,hy,nx,ny):
    Ml=np.zeros((nx*ny,nx*ny))
    for i in range(nx):
        for j in range(ny):
            im,ip,ii=IM[i],IP[i],II[i]
            jm,jp,jj=JM[j],JP[j],JJ[j]
            kl=(2/hx**2)*(np.cos(2*(np.pi/nx)*(row-1))-1)
            Ml[row, index(ii,jj)] = -2+kl*hy**2
            Ml[row, index(im,jj)] = 1
            Ml[row, index(ip,jj)] = 1
    Ml=(1/hy**2)*Ml



# en précalcul j'assemble le système linéaire M
#M=BuildPoisson(hx,hy,Nx,Ny)
# en précalcul on calcule l'inverse de M
#Minv=np.linalg.inv(M)
#print(np.matmul(Minv, M))

Ml=BuildTridiag(hx,yx,Nx,Ny)
Mlinv=np.linalg.inv(Ml)


nu=1e-3
dt = 0.25*hx**2 # stabilité
#
xc,yc= np.meshgrid(Xc,Yc, indexing='ij')
xu,yu= np.meshgrid(Xf,Yc, indexing='ij')
xv,yv= np.meshgrid(Xc,Yf, indexing='ij')

u,v=init_field()

t,Tmax=0,1.05

while(t<Tmax):
    t=t+dt
    NLu,NLv=get_NonLinear_terms(u,v)
    Lu,Lv=get_Linear_terms(u,v)
    dpdx,dpdy=get_gradient(pres)

    u_star = u + dt*( - dpdx - NLu + nu*Lu )
    v_star = v + dt*( - dpdy - NLv + nu*Lv )
    
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


    Gx,Gy=get_gradient(fi)
    u = u_star - dt*Gx
    v = v_star - dt*Gy
    pres = pres + fi
    div=get_divergence(u,v)

    print("%f %e"%(t,np.max(np.abs(div))))

omega=np.copy(u)
for i in II: 
    for j in JJ:
        omega[i,j]= (v[i,j]-v[IM[i],j])/hx-(u[i,j]-u[i,JM[j]])/hy

tp.tecplot_writer('test2d.dat', {'u': u,'v': u,'pres': omega}, Xf, Yf)

#matplotlib.pyplot.contourf(*args, data=None, **kwargs)
