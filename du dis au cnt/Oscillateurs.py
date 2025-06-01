
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv

T = 2*np.pi*100 # time of the simulation
T = 100.
M = 4000    # number of subdivisation
h = T/M
t = np.linspace(0,T,M) # definition of the array containing time-step

N = 32
q_old , p_old = np.zeros(N) , np.zeros(N) # poistion vitesse
q , p = np.zeros(N) , np.zeros(N)

q_animate = np.zeros((M,N))

f = open("Chaine_ES.dat","w")
# Euler symplectique
#q_old[0]=0
#p_old[0]=1.

#p_old=np.random.rand(N)
#p=np.random.rand(N)
p[0]=1
cpt=0

#q_hat=np.fft.fft(q)
#p[0]=1.

for tc in t:
    f.write("%f"%(tc))

    for n in np.arange(N):
        q[n] = q[n] + h*p[n]

#    p[0  ] = p[0  ] - h*( - q[N-1] + 2*q[0  ] - q[1] )
    p[0  ] = p[0  ] - h*(  + 2*q[0  ] - q[1] )
    for n in np.arange(1,N-1):
        p[n] = p[n] - h*( -q[n-1] + 2*q[n] - q[n+1])
#    p[N-1] = p[N-1] - h*( - q[N-2] + 2*q[N-1] - q[0] )
    p[N-1] = p[N-1] - h*( - q[N-2] + 2*q[N-1] )
    
    # calcul de l'Hamiltonien
    em = np.dot(p,p)
    em = em + q[0]*( -q[N-1] + 2*q[0] - q[1])
    em = em + q[N-1]*(-q[N-2] + 2*q[N-1] - q[0])
#    em = em + q[0]*(  + 2*q[0] - q[1])
#    em = em + q[N-1]*(-q[N-2] + 2*q[N-1])

    for n in np.arange(1,N-1):
        em = em + (-q[n-1] + 2*q[n] - q[n+1])*q[n]

    if (cpt==0):
        em0=em

    print( tc,np.abs(em-em0)  )
    for n in np.arange(N):
        f.write(" %e "%(p[n]))
    f.write("\n")
#    q_hat=np.fft.fft(q)
#    p_old=p
#    q_old=q
 #   q_animate[cpt,:]=q[:]
    cpt=cpt+1
    
f.close()

# animate the data


exit()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

k = 2*np.pi
w = 2*np.pi
dt = 0.01

xmin = 0
xmax = 3
nbx = 151

x = np.linspace(xmin, xmax, nbx)

fig = plt.figure() # initialise la figure
line, = plt.plot([], []) 
plt.xlim(xmin, xmax)
plt.ylim(-1, 1)

def init():
    line.set_data([], [])
    return line,

def animate(i): 
    t = i * dt
    y = np.cos(k*x - w*t)
    line.set_data(x, y)
    return line,
 
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=100, 
                              interval=1, blit=True, repeat=False)

plt.show()
exit()

