import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

L=4
lam=4
s=2
sp=2
h=5
delta_x=L/4
n=100
T_inf=0
Tmur=400

c=lam*s/delta_x


T1=np.zeros(n)
T2=np.zeros(n)
T3=np.zeros(n)
T4=np.zeros(n)


for i in range (n-1):
    T1[i+1]= c*(T2[i]-2*T1[i]+Tmur) -h*2*sp*(T1[i]-T_inf)
    T2[i+1]= c*(T3[i]-2*T2[i]+T1[i]) -h*2*sp*(T2[i]-T_inf)
    T3[i+1]= c*(T4[i]-2*T3[i]+T2[i]) -h*2*sp*(T3[i]-T_inf)
    T4[i+1]= c*(T4[i]-T3[i]) -h*(2*sp+s)*(T4[i]-T_inf)

plt.plot(n,T1)


