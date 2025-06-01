import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as s

#Exercice 2

#euler explicite

#Je prepare les pas
p={0:0.05,1:0.03,2:0.01}
#J'execute euler explicite avec les differents pas
for j in range (3):
    X=np.arange(0,2+p[j],p[j]) #vecteur des temps
    Y=np.zeros(len(X)) #future vecteur des solutions
    Y[0]=0 #ligne inutile dans ce cas
    for i in range(len(X)-1):
        Y[i+1]=Y[i]-50*p[j]*(Y[i]-m.cos(X[i]))
    plt.figure(j)
    plt.plot(X,Y)
    plt.show()
#La solution explose pour un pas grand mais converge pour un pas inférieur à 0.02. (convergence conditionelle)

#le meme test entre 0 et 100
X4=np.arange(0,100+0.01,0.01)
Y4=np.zeros(len(X4))
Y4[0]=0 
for i in range(len(X4)-1):
    Y4[i+1]=Y4[i]-50*0.01*(Y4[i]-m.cos(X4[i]))
plt.figure(4)
plt.plot(X4,Y4)
plt.show()


#euler implicite sans Newton

for j in range (3):
    X=np.arange(0,2+p[j],p[j]) #vecteur des temps
    Y=np.zeros(len(X)) #future vecteur des solutions
    Y[0]=0 #ligne inutile dans ce cas
    for i in range(len(X)-1):
        Y[i+1]=(Y[i]+50*p[j]*m.cos(X[i]))/51
    plt.figure(j+4)
    plt.plot(X,Y)
    plt.show()

#le meme test entre 0 et 100
X4=np.arange(0,100+0.01,0.01)
Y4=np.zeros(len(X4))
Y4[0]=0 
for i in range(len(X4)-1):
    Y4[i+1]=(Y4[i]+50*0.01*m.cos(X4[i]))/51
plt.figure(8)
plt.plot(X4,Y4)
plt.show()

#euler implicite avec Newton


for j in range (3):
    X=np.arange(0,2+p[j],p[j])
    Y=np.zeros(len(X))
    Y[0]=0 
    for i in range(len(X)-1):
        Y[i+1]=s.newton(lambda z:-50*(z-m.cos(X[i])),Y[i])
    plt.figure(j+8)
    plt.plot(X,Y)
    plt.show()

#Exercice 3

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
#couleurs
c={0:'blue',1:'orange',2:'red'} #respectivement S I et R

plt.figure(49)
for i in range (3):
    plt.plot(X,Y[:,i],color=c[i])
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
        plt.plot(X,Y[:,i],color=c[i])
    plt.show()
    print(nmbr)


#on voit que Sn+In+Rn=S0+I0+R0 donc la verification numerique est ok


#euler implicite

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
    Y[i+1,2]=s.newton(lambda z: a*z,Y[i,2],lambda zprime: zprime) #Je commence avec R cette fois
    Y[i+1,1]=s.newton(lambda z3: r*Y[i+1,2]*z3 - a*Y[i+1,2],Y[i,1],lambda zprime3: r*Y[i+1,2]) #I
    Y[i+1,0]=s.newton(lambda z2: -r*Y[i+1,2]*z2,Y[i,0],lambda zprime2: -r*Y[i+1,2] ) #S
    

plt.figure(49)
for i in range (3):
    plt.plot(X,Y[:,i],color=c[i])
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
        Y[i+1,2]=s.newton(lambda z: a*z,Y[i,2],lambda zprime: zprime) #Je commence avec R cette fois pour l'utiliser pour le calcul de S et I
        Y[i+1,1]=s.newton(lambda z3: r*Y[i+1,2]*z3 - a*Y[i+1,2],Y[i,1],lambda zprime3: r*Y[i+1,2]) #I
        Y[i+1,0]=s.newton(lambda z2: -r*Y[i+1,2]*z2,Y[i,0],lambda zprime2: -r*Y[i+1,2] ) #S
        #J'ai choisi d'utiliser la variable nombre pour compter le nombre de fois où la somme est differente de S0+I0+R0
        test=Y[i,0]+Y[i,1]+Y[i,2]
        if(test!=S0+I0+R0):
            nmbr=+1
    
    #ATTENTION beacoup de graphs
    plt.figure(i+50)
    for i in range (3):
        plt.plot(X,Y[:,i],color=c[i])
    plt.show()
    print(nmbr)


#on voit que Sn+In+Rn=S0+I0+R0 donc la verification numerique est ok