import numpy as np
import matplotlib.pyplot as plt
import math as m


#euler explicite
#pas superieur a 0.04
X=np.arange(0,2+0.05,0.05)
Y=np.zeros(len(X))
Y[0]=0 #ligne inutile
for i in range(len(X)-1):
    Y[i+1]=Y[i]-50*0.05*(Y[i]-m.cos(X[i]))
plt.figure(1)
plt.plot(X,Y)
plt.show()
    
#pas entre 0.02 et 0.04
X2=np.arange(0,2+0.03,0.03)
Y2=np.zeros(len(X2))
Y2[0]=0 #ligne inutile
for i in range(len(X)-1):
    Y2[i+1]=Y2[i]-50*0.03*(Y2[i]-m.cos(X2[i]))
plt.figure(2)
plt.plot(X2,Y2)
plt.show()

#pas inferieur a 0.02
X3=np.arange(0,2+0.01,0.01)
Y3=np.zeros(len(X3))
Y3[0]=0 #ligne inutile
for i in range(len(X3)-1):
    Y3[i+1]=Y3[i]-50*0.01*(Y3[i]-m.cos(X3[i]))
plt.figure(3)
plt.plot(X3,Y3)
plt.show()


#La solution explose pour un pas grand mais converge pour un pas inférieur à 0.02

#le même test entre 0 et 100
X4=np.arange(0,100+0.01,0.01)
Y4=np.zeros(len(X4))
Y4[0]=0 #ligne inutile
for i in range(len(X4)-1):
    Y4[i+1]=Y4[i]-50*0.01*(Y4[i]-m.cos(X4[i]))
plt.figure(4)
plt.plot(X4,Y4)
plt.show()

#on voit bien le caractere ossilant de la solution


#euler implicite
#pas superieur a 0.04
X=np.arange(0,2+0.05,0.05)
Y=np.zeros(len(X))
Y[0]=0 #ligne inutile
for i in range(len(X)-1):
    Y[i+1]=(Y[i]+50*0.05*m.cos(X[i]))/51
plt.figure(12)
plt.plot(X,Y)
plt.show()
    
#pas entre 0.02 et 0.04
X2=np.arange(0,2+0.03,0.03)
Y2=np.zeros(len(X2))
Y2[0]=0 #ligne inutile
for i in range(len(X)-1):
    Y2[i+1]=(Y2[i]+50*0.03*m.cos(X2[i]))/51
plt.figure(22)
plt.plot(X2,Y2)
plt.show()

#pas inferieur a 0.02
X3=np.arange(0,2+0.01,0.01)
Y3=np.zeros(len(X3))
Y3[0]=0 #ligne inutile
for i in range(len(X3)-1):
    Y3[i+1]=(Y3[i]+50*0.01*m.cos(X3[i]))/51
plt.figure(32)
plt.plot(X3,Y3)
plt.show()


#le même test entre 0 et 100
X4=np.arange(0,100+0.01,0.01)
Y4=np.zeros(len(X4))
Y4[0]=0 #ligne inutile
for i in range(len(X4)-1):
    Y4[i+1]=(Y4[i]+50*0.01*m.cos(X4[i]))/51
plt.figure(42)
plt.plot(X4,Y4)
plt.show()

