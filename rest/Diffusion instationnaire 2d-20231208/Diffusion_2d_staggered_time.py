import numpy as np
import matplotlib.pyplot as plt


# le maillage est défini ici
xmin , xmax , nx = 0. , 1. , 16
ymin , ymax , ny = 0. , 1. , 16


# des fonctions pour définir le terme source
def udf_s(x,y):
    k=2
    res= +np.sin(2*np.pi*k*x)*np.sin(2*np.pi*k*y)*(2*np.pi)**2*2
#    res= 0
    return res

# des fonctions pour définir les valeurs sur les bords
def udf_u(x,y):
    return 0.
def udf_gauche(x,y):
    return udf_u(x,y)
def udf_droit(x,y):
    return udf_u(x,y)
def udf_bas(x,y):
    return udf_u(x,y)
def udf_haut(x,y):
    return udf_u(x,y)


# le pas d'espace pour x et y
hx , hy = (xmax-xmin)/nx , (ymax-ymin)/ny
# les coordonnées en xc et yc du centre des cellules
xc = np.linspace(xmin+.5*hx,xmax-0.5*hx,nx)
yc = np.linspace(ymin+.5*hy,ymax-0.5*hy,ny)

# c'est juste pour éviter d'avoir à réécrire les indexes dans les boucles
ii,jj=np.arange(nx),np.arange(ny)

nxy=nx*ny # la dimension du vecteur des inconnues
M=np.zeros((nxy,nxy)) # le système linéaire
b=np.zeros((nxy)) # le second membre
x=np.zeros((nxy)) # la solution
I=np.eye((nxy)) # la solution

u,u_old=np.zeros((nxy)),np.zeros((nxy))
rhs=np.zeros((nxy)) 

def contribution_implicit(xc,yc):
    M=np.zeros((nxy,nxy)) # le système linéaire
    rhs=np.zeros(nxy) # le système linéaire
    for j in jj:
        for i in ii:
            row = i + j*(nx) # ligne de la matrice à remplir
            
            Ap =  i + j*(nx) # cellule courante rem Ap == row :-)
            Aw,Ae = i-1 + j*nx , i+1 + j*nx #
            M[row,Ap] = M[row,Ap] - 2/hx**2 # j'ajoute la contribution diagonale
            if (i==0) : # je suis sur le bord gauche
                M[row,Ap] =  M[row,Ap] - 1/hx**2
                M[row,Ae] =  M[row,Ae] + 1/hx**2
                rhs[row] =  rhs[row] + 2*udf_gauche( xmin , yc[j]  )/hx**2
            elif (i==(nx-1)): # je suis sur le bord droit
                M[row,Ap] =  M[row,Ap] - 1/hx**2
                M[row,Aw] =  M[row,Aw] + 1/hx**2
                rhs[row] =  rhs[row] + 2*udf_droit( xmax , yc[j]  )/hx**2
            else: # cellule intérieure en pour ox
                M[row,Aw] =  M[row,Aw] + 1/hx**2
                M[row,Ae] =  M[row,Ae] + 1/hx**2

                # contribution dans la direction oy
            As,An = i + (j-1)*nx , i + (j+1)*nx

            M[row,Ap] = M[row,Ap] - 2/hy**2  # j'ajoute la contribution diagonale
            if (j==0) : # je suis sur le bord bas
                M[row,Ap] =  M[row,Ap] - 1/hy**2
                M[row,An] =  M[row,An] + 1/hy**2
                rhs[row] =  rhs[row] + 2*udf_bas( xc[i] , ymin  )/hy**2

            elif (j==(ny-1)): # je suis sur le bord haut
                M[row,Ap] =  M[row,Ap] - 1/hy**2
                M[row,As] =  M[row,As] + 1/hy**2
                rhs[row] =  rhs[row] + 2*udf_haut( xc[i] , ymax  )/hy**2
                
            else:  # cellule intérieur en pour oy
                M[row,As] =  M[row,As] + 1/hy**2
                M[row,An] =  M[row,An] + 1/hy**2

    return (M,rhs)

def contribution_explicit(xc,yc,u_old):
    rhs_exp=np.zeros(nxy) # le système linéaire
    for j in jj:
        for i in ii:
            row = i + j*(nx) # ligne de la matrice à remplir
            Ap =  i + j*(nx) # cellule courante rem Ap == row :-)

            Aw,Ae = i-1 + j*nx , i+1 + j*nx #
            rhs_exp[row] =  -2/hx**2*u_old[Ap] # j'ajoute la contribution diagonale
            if (i==0) : # je suis sur le bord gauche
                rhs_exp[row] = rhs_exp[row] - 1/hx**2*u_old[Ap]
                rhs_exp[row] = rhs_exp[row] + 1/hx**2*u_old[Ae]
                rhs_exp[row] = rhs_exp[row] + 2*udf_gauche( xmin , yc[j]  )/hx**2
            elif (i==(nx-1)): # je suis sur le bord droit
                rhs_exp[row] = rhs_exp[row] - 1/hx**2*u_old[Ap]
                rhs_exp[row] = rhs_exp[row] + 1/hx**2*u_old[Aw]
                rhs_exp[row] = rhs_exp[row] + 2*udf_droit( xmax , yc[j]  )/hx**2
            else: # cellule intérieur en pour ox
                rhs_exp[row] = rhs_exp[row] + 1/hx**2*u_old[Aw]
                rhs_exp[row] = rhs_exp[row] + 1/hx**2*u_old[Ae]

            As,An = i + (j-1)*nx , i + (j+1)*nx #
            rhs_exp[row] = rhs_exp[row] - 2/hy**2*u_old[Ap] # j'ajoute la contribution diagonale
            if (j==0) : # je suis sur le bord bas
                rhs_exp[row] = rhs_exp[row] - 1/hy**2*u_old[Ap]
                rhs_exp[row] = rhs_exp[row] + 1/hy**2*u_old[An]
                rhs_exp[row] = rhs_exp[row] + 2*udf_bas( xc[i] , ymin  )/hy**2
            elif (j==(ny-1)): # je suis sur le bord haut
                rhs_exp[row] = rhs_exp[row] - 1/hy**2*u_old[Ap]
                rhs_exp[row] = rhs_exp[row] + 1/hy**2*u_old[As]
                rhs_exp[row] = rhs_exp[row] + 2*udf_haut( yc[i] , ymax )/hy**2
            else: # cellule intérieur en pour oy
                rhs_exp[row] = rhs_exp[row] + 1/hy**2*u_old[As]
                rhs_exp[row] = rhs_exp[row] + 1/hy**2*u_old[An]

    return rhs_exp



# calcul symbolique ---------------------
import sympy as sym
xx = sym.Symbol('x') # xx est un symbole
yy = sym.Symbol('y') # a est un symbole
tt = sym.Symbol('t') # a est un symbole

u_mms = xx+yy # on construit une expression analytique
e=1
# puis on calcule le Laplacien
dg = sym.diff(u_mms,tt,1) - (sym.diff(u_mms,xx,2) + sym.diff(u_mms,yy,2))
print(dg)
# et on transforme en fonction numérique
scm_mms=sym.lambdify([tt,xx,yy],dg)
u_mms=sym.lambdify([tt,xx,yy],u_mms)
#---------------------------------------
exit()

#
theta = 0.5
dt = 1e-1
Tmax = 10*dt
tc = 0
nu = 1
while(tc<Tmax):
    tc= tc+dt
    # on initialise le terme source
    for j in jj:
        for i in ii:
            rhs[i+j*nx] = udf_s( xc[i] , yc[j] ) #
    # je contruis le système linéaire
    M,rhs_im = contribution_implicit(xc,yc) # on assemble le laplacien
    M = I - dt*nu*M * theta # on discrétise en temps 
    rhs_im = -dt*nu*rhs_im * theta 

    rhs_ex = contribution_explicit(xc,yc,u_old) # on cacule les terme explicites
    rhs_ex = - dt*nu*rhs_ex * (1-theta)
    
#    b = u_old - rhs_im - rhs_ex  + rhs
    b = u_old + rhs*dt - rhs_im - rhs_ex
    u = np.linalg.solve(M,b)
    u_old = u

    print(np.min(u),np.max(u))







# pour ploter 
fi=np.zeros((nx,ny))
for j in jj:
    for i in ii:
        row = i + j*(nx)
        #fi[i,j] = np.abs(u[row]-udf_u(xc[i],yc[j]))
        fi[i,j] = u[row]
print(np.max(fi))

print("----")
for j in jj:
    aa=udf_u(xc[3],yc[j])
    print(("%f %f ")%(u[3 + j*nx]-aa,u[3 + j*nx]-aa))



        
xc_v, yc_v = np.meshgrid(xc, yc,indexing='ij')
plt.contourf(xc_v,yc_v,fi)
plt.show()

for j in jj:
    for i in ii:
        row = i + j*(nx)
        fi[i,j] = rhs_im[row]

