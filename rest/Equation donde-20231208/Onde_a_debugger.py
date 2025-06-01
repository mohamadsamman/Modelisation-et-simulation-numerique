import numpy as np
import matplotlib.pyplot as plt


r_min , r_max  = 0.0 , 1. # domaine [r_min,r_max]x[0,2pi]
Nr , Na = 64 , 64 # nd de cellules en r et \theta

hr , ha = (r_max-r_min)/Nr , 2*np.pi/Na # taille des cellules en 
r = np.linspace(r_min,r_max,Nr) # r_i={}
a = np.linspace(0,2*np.pi,Na,endpoint=False) # theta_j={}

ii , jj = np.arange(Nr),np.arange(Na) # espace des indices en r et \theta

print(a)


Lu=np.zeros((Nr,Na))
jm,jp=jj-1,jj+1
jm[0],jp[-1]=Na-1,0
print(jp)

def compte_rhs(u):
    Lu=np.zeros((Nr,Na))
    for j in np.arange(0,Na):

        for i in np.arange(1,Nr-1):
            Lu[i,j]=(u[i+1,j]-2*u[i,j]+u[i-1,j])/hr**2 + r[i]**(-1)*(u[i+1,j]-u[i-1,j])/(2*hr) \
              + r[i]**(-2)*(u[i,jp[j]]-2*u[i,j]+u[i,jm[j]])/ha**2
        # conditions aux limites de i=Nr-1
        i=Nr-1
        Lu[i,j]=(u[i-1,j]-2*u[i,j]+u[i-1,j])/hr**2 + r[i]**(-1)*(u[i-1,j]-u[i-1,j])/(2*hr) \
          + r[i]**(-2)*(u[i,jp[j]]-2*u[i,j]+u[i,jm[j]])/ha**2
        i=0
        u0_mean=sum(u[1,:])/Na
        print(u0_mean)
        # conditions aux limites de i=0
        Lu[i,j]=4*(u0_mean-u[0,j])/hr**2
    return Lu

un=np.zeros((Nr,Na))
u =np.zeros((Nr,Na))
uo=np.zeros((Nr,Na))

dt=1e-3
tmax=1000*dt
tc=0
#u[0,:]=1
u[32,1]=1
uo=u
while (tc<tmax):
    tc = tc + dt
    Lu = compte_rhs(u)
    un = ( 2*u - uo ) + dt**2*Lu
    uo = u
    u = un
    
rv, av = np.meshgrid(r,a,indexing='ij')
x,y=rv*np.cos(av),rv*np.sin(av)
plt.contour(x,y,u)
plt.show()
