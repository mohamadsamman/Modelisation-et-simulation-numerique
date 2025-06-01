import numpy as np
import matplotlib.pyplot as plt


ri,ro=1.,2. # rayon interieur, extérieur
Nr,Na = 32,64    # nombre d'éléments
dr,da = (ro-ri)/(Nr+1) , 2*np.pi/Na  

#   ri          ro  
#   o--x--x--x--o   Nr= 3 noeuds intérieurs + 2 frontiéres == 4 cellules
rf = np.linspace(ri+dr,ro-dr,Nr)

#   0           2pi  
#   x--x--x--x--o   Na= 4 noeuds intérieurs 
af = np.linspace(0,2*np.pi-da,Na)

N=Nr*Na
M=np.zeros((N,N))
scm=-np.ones(N)

# calcul symbolique ---------------------
import sympy as sym
rr = sym.Symbol('rr') # r est un symbole
aa = sym.Symbol('aa') # a est un symbole
RR = sym.Symbol('RR') # a est un symbole

u_mms = -(rr-ri)*(rr-ro) # on construit une expression analytique
RR = (rr-ri)*(ro-ri) # on construit une expression analytique

u_mms=sym.sin(2*np.pi*RR)*sym.cos(4*aa)
e=1
# puis on calcule le Laplacien
dg = sym.diff(u_mms,rr,2) + 1/rr * sym.diff(u_mms,rr)*e + 1/rr**2 * sym.diff(u_mms,aa,2)
print(dg)
# et on transforme en fonction numérique 
scm_mms=sym.lambdify([rr,aa],dg)
u_mms=sym.lambdify([rr,aa],u_mms)
#---------------------------------------

for j in np.arange(Na):
    for i in np.arange(Nr):
        p=i+j*Nr
        scm[p] = scm_mms( rf[i] , af[j] )


# itération sur les noeuds intérieurs
for j in np.arange(Na):
    for i in np.arange(Nr):
        
        p = i + j*Nr
        if (i==0):
            idx = [p,p+1]
            ste = [-2/dr**2 , 1/dr**2 + 1/(2*dr*rf[i])*e]
        elif (i==Nr-1):
            idx = [ p-1 , p ]
            ste = [ 1/dr**2 - 1/(2*dr*rf[i])*e , -2/dr**2 ]
        else:
            idx = [p-1,p,p+1]
            ste = [ 1/dr**2 - 1/(2*dr*rf[i])*e , -2/dr**2 , 1/dr**2 +1/(2*dr*rf[i])*e ]

        M[p, idx] += ste
        
        if (j==0):
            As,An = i + (Na-1)*Nr,i + (j+1)*Nr
            idx = [As,p,An]
            ste = [+1/da**2/rf[i]**2 , -2/da**2/rf[i]**2 , +1/da**2/rf[i]**2 ]
        elif (j==Na-1):
            As,An = i + (j-1)*Nr,i + 0*Nr
            idx = [As,p,An]
            ste = [+1/da**2/rf[i]**2 , -2/da**2/rf[i]**2 , +1/da**2/rf[i]**2 ]
        else:
            As,An = i + (j-1)*Nr,i + (j+1)*Nr
            idx = [As,p,An]
            ste = [+1/da**2/rf[i]**2 , -2/da**2/rf[i]**2 , +1/da**2/rf[i]**2 ]
        M[p, idx] += ste

U=np.linalg.solve(M,scm)



# du vecteur vers un champ 2d
fi=np.zeros((Nr,Na))
err=np.zeros((Nr,Na))
for j in np.arange(Na):
    for i in np.arange(Nr): 
        err[i,j]=U[i+j*Nr]-u_mms(rf[i],af[j])
        fi[i,j]=U[i+j*Nr]

print("err %f :"%np.max(np.abs(err)))
r, a = np.meshgrid(rf,af,indexing='ij')
x,y=r*np.cos(a),r*np.sin(a)

plt.contourf(x,y,fi)
plt.show()
exit()
