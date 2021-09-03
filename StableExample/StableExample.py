# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 14:50:16 2021

@author: bohne
"""


import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'ieee'])

T = np.linspace(-2000,2000, 50000)

H0 = 5.5E-6
alpha = 1/200
Mpl = 1

#More adjustable parameters
cF = 1
c30 = 0

def H(t):
    return H0*np.arctan(alpha*t)

def dH(t):
    return H0*alpha/(1+ alpha**2*t**2)
    
def G2F(x, c20,c21,c22,c31):
    return cF

def G2(x, c20,c21,c22,c31):
    return c20 + c21*x+c22*x**2

def dG2(x, c20,c21,c22,c31):
    return c21 + 2*c22*x

def ddG2(x, c20,c21,c22,c31):
    return 2*c22
    
def G3(x, c20,c21,c22,c31):
    return c30 + c31*x
    
def dG3(x, c20,c21,c22,c31):
    return c31

def ddG3(x, c20,c21,c22,c31):
    return 0

def G4(x, c20,c21,c22,c31):
    return Mpl**2/2
    
def A0p(t, c20,c21,c22,c31):
    return (3*c31*H(t)+np.sqrt(9*c31**2*H(t)**2 - 4*c21*c22))/(2*c22)

def A0m(t, c20,c21,c22,c31):
    return (3*c31*H(t)-np.sqrt(9*c31**2*H(t)**2 - 4*c21*c22))/(2*c22)

def dA0p(t, c20,c21,c22,c31):
    return 3*c31/(2*c22)*(1+ 3*c31*H(t)/(np.sqrt(9*c31**2*H(t)**2 - 4*c21*c22)))*dH(t)

def dA0m(t, c20,c21,c22,c31):
    return 3*c31/(2*c22)*(1-3*c31*H(t)/(np.sqrt(9*c31**2*H(t)**2 - 4*c21*c22)))*dH(t)

##Coefficients
def qt(t, A0, c20,c21,c22,c31):
    x = A0**2/2
    return 2*G4(x, c20,c21,c22,c31)

def Ft(t, A0, c20,c21,c22,c31):
    x = A0**2/2
    return 2*G4(x, c20,c21,c22,c31)

def qv(t, A0, c20,c21,c22,c31):
    x = A0**2/2
    return G2F(x, c20,c21,c22,c31)

def Fv(t, A0, c20,c21,c22,c31):
    x = A0**2/2
    return G2F(x, c20,c21,c22,c31)

def w1(t,A0, c20,c21,c22,c31):
    x = A0**2/2
    return -4*G4(x, c20,c21,c22,c31)*H(t)+A0**3*dG3(x, c20,c21,c22,c31)


def w2(t,A0, c20,c21,c22,c31):
    return w1(t, A0, c20,c21,c22,c31) + 2*H(t)*qt(t, A0, c20,c21,c22,c31)

def w3(t, A0, c20,c21,c22,c31):
    return -2*A0**2*qv(t, A0, c20,c21,c22,c31)

def w4(t, A0, c20,c21,c22,c31):
    x = A0**2/2
    return 1/2*A0**4*ddG2(x, c20,c21,c22,c31)+3/2*A0**3*H(t)*dG3(x, c20,c21,c22,c31) \
        -3/2*A0**5*H(t)*ddG3(x, c20,c21,c22,c31)-6*G4(x, c20,c21,c22,c31)*H(t)**2

def w5(t, A0, c20,c21,c22,c31):
    return w4(t, A0, c20,c21,c22,c31)-3/2*H(t)*(w1(t,A0, c20,c21,c22,c31)+w2(t,A0, c20,c21,c22,c31))

def w6(t,A0, c20,c21,c22,c31):
    x = A0**2/2
    return -A0**2*dG3(x, c20,c21,c22,c31)

def w7(t,A0,dA0, c20,c21,c22,c31):
    x = A0**2/2
    return dA0*dG3(x, c20,c21,c22,c31)

def w8(t,A0, c20,c21,c22,c31):
    return 3*H(t)*w1(t, A0, c20,c21,c22,c31) - 2*w4(t, A0, c20,c21,c22,c31)

def qs(t, A0, c20,c21,c22,c31):
    return qt(t,A0, c20,c21,c22,c31)*(3*w2(t,A0, c20,c21,c22,c31)**2+4*w5(t,A0, c20,c21,c22,c31)*qt(t,A0, c20,c21,c22,c31))/(2*H(t)*qt(t,A0, c20,c21,c22,c31) +w2(t,A0, c20,c21,c22,c31))**2


def Fs(t,A0, dA0,c20t,c21t,c22t,c31t):
    x = A0**2/2
    return -2*A0**2*G4(x,c20t,c21t,c22t,c31t)*(4*A0*G4(x,c20t,c21t,c22t,c31t)*H(t)*dG3(x,c20t,c21t,c22t,c31t)*G2F(x,c20t,c21t,c22t,c31t) \
        +8*G4(x,c20t,c21t,c22t,c31t)*dA0*dG3(x,c20t,c21t,c22t,c31t)*G2F(x,c20t,c21t,c22t,c31t) \
            +A0**4*dG3(x,c20t,c21t,c22t,c31t)**2*G2F(x,c20t,c21t,c22t,c31t)\
    - 4*A0**2*G4(x,c20t,c21t,c22t,c31t)*(dG3(x,c20t,c21t,c22t,c31t)**2 \
    - dA0*ddG3(x,c20t,c21t,c22t,c31t)*G2F(x,c20t,c21t,c22t,c31t)))/(G2F(x,c20t,c21t,c22t,c31t)*(4*G4(x,c20t,c21t,c22t,c31t)*H(t) \
                + A0**3*dG3(x,c20t,c21t,c22t,c31t)))**2
        
def Rho(t,A0,c20t,c21t,c22t,c31t):
    x = A0**2/2
    return G2(x,c20t,c21t,c22t,c31t) + 6*G4(x,c20t,c21t,c22t,c31t)*H(t)**2

def p(t,A0, dA0,c20t,c21t,c22t,c31t):
    x = A0**2/2
    return -G2(x,c20t,c21t,c22t,c31t)-A0**2*dA0*dG3(x,c20t,c21t,c22t,c31t)-6*G4(x,c20t,c21t,c22t,c31t)*H(t)**2\
        -4*G4(x,c20t,c21t,c22t,c31t)*dH(t)
 

c20 = 0
c21 = 0.75
c22 = -1
c31 = 3       
 
fig,ax=plt.subplots()

plt.plot(T, H(T), color = "green", linestyle = "--", label = "$\ H(t)$")
a = '\\cdot\dot{H}(t)'
ax.plot(T, 300*dH(T), color = "orange", linestyle = ":",label = '300$%s$'%a)
ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.legend()
plt.savefig("plots/Fig51_H.png")

fig,ax=plt.subplots()
plt.plot(T, -A0p(T, c20,c21,c22,c31), color = "black", linestyle = "-", label = "$\ -A_0(t)^{+}$")
ax.plot(T, A0m(T, c20,c21,c22,c31), color = "red", linestyle = "--", label = "$\ A_0(t)^{-}$")
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(0.865916, 0.866135)
plt.legend()
plt.savefig("plots/Fig52_A0.png")

A0pArray = A0p(T, c20,c21,c22,c31)
A0mArray = A0m(T, c20,c21,c22,c31)


Rhop = np.zeros(len(T))
Rhom = np.zeros(len(T))
for i,t in enumerate(T):
    Atemp = A0p(t, c20,c21,c22,c31)
    Rhop[i] = Rho(t,Atemp, c20,c21,c22,c31)
    
    Atemp = A0m(t, c20,c21,c22,c31)
    Rhom[i] = Rho(t,Atemp, c20,c21,c22,c31)
    
RhoMin = min(Rhop)
RhoMax = max(Rhop)
RhoDiff = RhoMax- RhoMin

fig,ax=plt.subplots()
plt.plot(T, Rhop, label = "$\\rho_M^{+}$")
plt.plot(T, Rhom, label = "$\\rho_M^{-}$")
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(RhoMin- 0.2*RhoDiff,RhoMax+0.3*RhoDiff )
plt.legend()
plt.savefig("plots/Fig53c_rho.png")


RhoPp = np.zeros(len(T))
RhoPm = np.zeros(len(T))
for i,t in enumerate(T):
    Atemp = A0p(t, c20,c21,c22,c31)
    dAtemp = dA0p(t, c20,c21,c22,c31)
    rhotemp = Rho(t,Atemp, c20,c21,c22,c31)
    ptemp = p(t,Atemp, dAtemp, c20,c21,c22,c31)
    RhoPp[i] = rhotemp+ptemp
    
    Atemp = A0m(t, c20,c21,c22,c31)
    dAtemp = dA0m(t, c20,c21,c22,c31)
    rhotemp = Rho(t,Atemp, c20,c21,c22,c31)
    ptemp = p(t,Atemp, dAtemp,c20,c21,c22,c31)
    RhoPm[i] = rhotemp+ptemp

fig,ax=plt.subplots()
plt.plot(T, Rhop, label = "$\\rho_M^{+}$")
plt.plot(T, Rhom, label = "$\\rho_M^{-}$")
plt.plot(T, RhoPp*1000000, label = "$\\rho_M^{+} + P_M^{+}$")
plt.plot(T, RhoPm*1000000, label = "$\\rho_M^{-} + P_M^{-}$")
ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.legend()


RhopPMax = max(RhoPp)
RhopPMin = min(RhoPp)
RhopPDiff = RhopPMax- RhopPMin

fig,ax=plt.subplots()
plt.plot(T, RhoPp, label = "$\\rho_M^{+} + P_M^{+}$")
plt.plot(T, RhoPm, label = "$\\rho_M^{-} + P_M^{-}$")
ax.axhline(y=0, color='k', linewidth = 0.5)
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(RhopPMin-0.2*RhopPDiff, RhopPMax + 0.3*RhopPDiff)
plt.legend()
plt.savefig("plots/Fig53d_rhop.png")


qsp = np.zeros(len(T))
qsm = np.zeros(len(T))

for i,t in enumerate(T):
    Atemp = A0p(t,c20,c21,c22,c31)
    qsp[i] = qs(t, Atemp,c20,c21,c22,c31)
    
    Atemp = A0m(t,c20,c21,c22,c31)
    qsm[i] = qs(t, Atemp,c20,c21,c22,c31)
qsMax = max(qsp)
qsMin = min(qsp)
qsDiff = qsMax-qsMin
fig,ax=plt.subplots()
plt.plot(T, qsp, label = "$\ q_S^{+}$")
plt.plot(T, qsm, label = "$\ q_S^{-}$")
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(qsMin-qsDiff, qsMax+qsDiff)
plt.legend()
plt.savefig("plots/Fig53a_qs.png")

Fsp = np.zeros(len(T))
Fsm = np.zeros(len(T))
for i,t in enumerate(T):
    Atemp = A0p(t,c20,c21,c22,c31)
    dAtemp = dA0p(t,c20,c21,c22,c31)
    Fsp[i] = Fs(t, Atemp, dAtemp,c20,c21,c22,c31)
    
    Atemp = A0m(t,c20,c21,c22,c31)
    dAtemp = dA0m(t,c20,c21,c22,c31)
    Fsm[i] = Fs(t, Atemp, dAtemp,c20,c21,c22,c31)

fig,ax=plt.subplots()
plt.plot(T, Fsp, label = "$\ F_S^{+}$")
plt.plot(T, Fsm, label = "$\ F_S^{-}$")
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(1.66611, 1.66722)
plt.legend()
plt.savefig("plots/Fs.png")


csp = np.zeros(len(T))
csm = np.zeros(len(T))
for i,t in enumerate(T):
    csp[i] = Fsp[i]/qsp[i]
    csm[i] = Fsm[i]/qsm[i]

fig,ax=plt.subplots()
plt.plot(T, csp, label = "$\ c_S^{2,+}$")
plt.plot(T, csm, label = "$\ c_S^{2,-}$")
ax.axvline(x=0, color='k', linewidth = 0.5)
plt.xlabel("$\ t$")
plt.ylim(0.692641, 0.691975)
plt.legend()
plt.savefig("plots/Fig53b_cs.png")


