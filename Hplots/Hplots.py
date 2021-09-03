# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 20:03:09 2021

@author: bohne
"""


import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'ieee'])

T = np.linspace(-100,1000, 50000)
d1 = 5.5*10**(-6)

def H(t):
    return d1*t/(1+d1*t**2)

def dH(t):
    return -2*d1**2*t**2/(1+d1*t**2)**2 + d1/(1 + d1*t**2)





def a(t):
    return np.sqrt(1 + d1*t**2)

def da(t):
    return d1*t/(np.sqrt(1+d1*t**2))


T = np.linspace(-10000,10000, 50000)
plt.figure()
plt.plot(T,500*a(T), label = "500")
plt.plot(T,250*a(T), label = "250")
plt.plot(T,750*a(T), label = "750")
plt.legend()


T = np.linspace(-2500,2500, 50000)



tau = 10
def H(t):
    return t/(3*(tau**2 + t**2))

def a(t):
    return (tau**2 + t**2)**(1/6)
T = np.linspace(-70,70, 50000)

def dH(t):
    return (tau**2- t**2)/(3*(tau**2+t**2)**2)


plt.figure()
plt.plot(T,H(T))

fig,ax=plt.subplots()

horizon3 = np.zeros(50000)
for i,t in enumerate(T):
    horizon3[i] = 1/abs(H(t))

ax.axvspan(-tau,tau,facecolor = 'magenta',alpha = 0.1, label ="Bounce", zorder = 1)
ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.plot(T, horizon3, color = "green", linestyle = "--", label = "$\ H(t)^{-1}$", zorder = 3)
plt.plot(T,40*a(T), label = '$\\lambda = a(t)/k$', zorder = 4)
plt.ylim([0,200])
ax.axvline(x = tau, color='magenta', linewidth = 0.7, zorder = 1)
ax.axvline(x = -tau, color='magenta', linewidth = 0.7, zorder = 1)
ax.axhline(y = 10, color = 'red', linewidth = 0.7, zorder = 2)
ax.axhspan(0,10,facecolor = 'red',alpha = 0.25, label ="Trans-Planckian", zorder = 2)

ax.set_yticklabels([])

plt.xlabel("$\ t$")
plt.legend(loc = (0.01,0.1))
plt.savefig("plots/Fig11_Horizon.png")

T = np.linspace(-45,45, 50000)
fig,ax=plt.subplots()
plt.plot(T, H(T), color = "green", linestyle = "--", label = "$\ H(t)$", zorder = 2)
a = '\\cdot\dot{H}(t)'
ax.axvline(x = tau, color='magenta', linewidth = 0.7, zorder = 1)
ax.axvline(x = -tau, color='magenta', linewidth = 0.7, zorder = 1)

ax.axvspan(-tau,tau,facecolor = 'magenta',alpha = 0.1, label ="Bounce", zorder = 1)
ax.plot(T, 5*dH(T), color = "orange", linestyle = ":",label = '5$%s$'%a, zorder = 3)


ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)


ax.set_xlabel("$\ t$")
ax.xaxis.set_label_coords(0.9,-0.1)

plt.legend()
plt.savefig("plots/Fig21_NEC_H.png")



T = np.linspace(-5,35, 50000)
fig,ax=plt.subplots()
plt.plot(T, H(T), color = "green", linestyle = "--", label = "$\ H(t)$", zorder = 2)
ax.plot(T, 5*dH(T), color = "orange", linestyle = ":",label = '5$%s$'%a, zorder = 3)
a = '\\cdot\dot{H}(t)'
ax.axvline(x = tau, color='magenta', linewidth = 0.7, zorder = 1, label = "$\ t^{\\star}$")

ax.axvspan(-5,tau,facecolor = 'magenta',alpha = 0.1, label ="Bounce", zorder = 1)



ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)

plt.xlim([-5,35])

ax.set_xlabel("$\ t$")

#plt.xlabel("$\ t$", loc = "right")
plt.legend(loc = "lower right")
plt.savefig("plots/Fig61_nogo_H.png")



