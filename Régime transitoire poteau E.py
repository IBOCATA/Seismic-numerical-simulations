# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 15:41:25 2024

@author: ilias
"""

import numpy as np
from scipy import stats
from scipy.stats import lognorm
from pylab import *
import matplotlib.pyplot as plt



# Message de début

print("********** START **********")
print("********** IBO ***********")
print("Projet generation des spectres de réponses")


# Number of samplepoints
N = 850
# sample spacing
T = 20 #durée enregistrée de séisme en s
temps = np.linspace(0.0, T, N)
#Données statistiques des séismes
#onde P
sigmap=0.6
mup=0
#onde S
sigmas=0.8
mus=0.05
#génération du signal sismique
d = np.concatenate((stats.norm.rvs(mup,sigmap,255),stats.norm.rvs(mus,sigmas,279),stats.norm.rvs(0.04,2.6,146),(stats.norm.rvs(mup,sigmap+0.1,170))))
# Tracé du signal sismique sample N
d=exp(-(temps-10)**2/10**2)*d


#Génération du spectre de réponse
#amortissement 
x=0.05
#incrément de temps
delta_T=temps[2]-temps[1]
plt.plot(temps, d)
plt.show() # Display plot to screen


#caractéristique de poteau en béton
E=21000000000
I=0.00045
L=3.5
M=5250
u=np.zeros((N))
v=np.zeros((N))
a=np.zeros((N))
w=sqrt(E*I/(L*M))
gamma=0.5
beta=0.25
for j in range(0,N-1):
    a[j+1]=(d[j+1]-w**2*u[j]-(2*x*w+2*delta_T*w**2)*v[j]-w*delta_T*(2*x*(1-gamma)+w*delta_T*(1-2*beta)/2)*a[j])/(1+((w*delta_T)**2)*beta+2*x*w*gamma)
    v[j+1]=v[j]+delta_T*((1-gamma)*a[j]+gamma*a[j+1])
    u[j+1]=u[j]+delta_T*v[j]+delta_T**2/2*((1-2*beta)*a[j]+2*beta*a[j+1])
time='temps en s'
plt.plot(temps,a)
plt.xlabel(time)
plt.ylabel('acceleration en m/s2')
plt.show()
plt.plot(temps,u)
plt.xlabel(time)
plt.ylabel('deplacement en m')
plt.show()
plt.plot(temps,v)
plt.xlabel(time)
plt.ylabel('vitesse en m/s')
plt.show()