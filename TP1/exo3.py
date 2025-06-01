import numpy as np
import matplotlib.pyplot as plt
import math as m


#initialisation des valeurs connues
S0=762
I0=1
R0=0
r=0.00218
a=0.44036
T=14

#euler explicite

#vecteur temps
X=np.linspace(0,14,1400)
#matrice de valeurs S,I,R (respectivement contenues dans les colones 0,1,2)
Y=np.zeros((len(X),3))
#initialistion du vecteur solution
Y[0,0]=S0
Y[0,1]=I0
Y[0,2]=R0
pas=X[3]-X[2]
for i in range(len(X)-1):
    Y[i+1,0]=Y[i,0]+pas*(-r*Y[i,0]*Y[i,1])
    Y[i+1,1]=Y[i,1]+pas*(r*Y[i,0]*Y[i,1]-Y[i,1])
    Y[i+1,2]=Y[i,2]+pas*(a*Y[i,1])


for i in range (3):
    plt.plot(X,Y[:,i])
plt.show()

#2
n=[14,20,30,40,50,60,70,80,90,100]

for i in range(len(n)):
    X=np.linspace(0,14,n[i])
    #matrice de valeurs S,I,R (respectivement contenues dans les colones 1,2,3)
    Y=np.zeros((len(X),3))
    #initialistion 
    Y[0,0]=S0
    Y[0,1]=I0
    Y[0,2]=R0
    nmbr=0
    for i in range(len(X)-1):
        Y[i+1,0]=Y[i,0]+pas*(-r*Y[i,0]*Y[i,1])
        Y[i+1,1]=Y[i,1]+pas*(r*Y[i,0]*Y[i,1]-Y[i,1])
        Y[i+1,2]=Y[i,2]+pas*(a*Y[i,1])
        #J'ai choisi d'utiliser la variable nombre pour compter le nombre de fois où la somme est differente de S0+I0+R0
        test=Y[i,0]+Y[i,1]+Y[i,2]
        if(test!=S0+I0+R0):
            nmbr=+1
    
    #ATTENTION beacoup de graphs
    plt.figure(i+50)
    for i in range (3):
        plt.plot(X,Y[:,i])
    plt.show()
    print(nmbr)


#on voit que Sn+In+Rn=S0+I0+R0 donc la verification numerique est ok
#La verification theorique vient du fait que S'+I'+R' = -rSI+rSI-aI+aI = 0 c'est à dire que la variation du nombre total d'individus est nulle. 