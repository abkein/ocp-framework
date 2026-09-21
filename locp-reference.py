from numba import njit, prange
from math import sqrt,pi
import numpy as np
import math
import random
import os
import time


m = 42500  # number of particles
mfv = 0    # number of core particles

Ndf = 150

G = 15.0
alpha = 0.0
fr = 0.09  # *sqrt(G/500)
tau0 = 0.01
beta = 1.9770761495708644
Cpe = 3.626300351593979e-11
wtime = 47.0
nch = 100

u0 = sqrt(2/G)
gamma = fr*tau0/2
delta = sqrt(2/G*fr/tau0)
Rb = (1 - alpha)*m**(1/3) - tau0*sqrt(1/8/pi/G)
R0 = m**(1/3)
Rcore = R0 - 3/2*2**(2/3)*(1/3)**0.5*(4*pi/3)**(1/3)
N = int(wtime/Cpe/m**beta)
N = 100*N + 1
if (N > 5400000):
    N = 5400001
# N *= 2 # no equilibrization

Na = int(abs(10000000*np.random.normal(0.0,1.0)))
if Na > 20*20000000:
    Na = 20000000
sn = 0.0
for i in range(Na):
    sn += np.sin(3.1415*i)

ini_file_name = f'ini_{m}_{G}.dat'
gini = np.zeros(0)
if os.path.exists(ini_file_name):
    idnum = 2
    os.remove(ini_file_name)
else:
    idnum = 1
    np.savetxt(ini_file_name, gini)

print()

f0 = np.zeros((m,3))   # total force on i-th particle at step 0, its j-th component
f1 = np.zeros((m,3))   # total force on i-th particle at step 1, its j-th component
fc0 = np.zeros((m,3))  # Coulomb force on i-th particle at step 0, its j-th component
fc1 = np.zeros((m,3))  # Coulomb force on i-th particle at step 1, its j-th component
fnd0 = np.zeros((m,3)) # total force on i-th particle at time step 0
fnd1 = np.zeros((m,3)) # total force on i-th particle at time step 1
x0 = np.zeros((m,3))   # j-th coordinate of i-th particle at step 0
x1 = np.zeros((m,3))   # j-th coordinate of i-th particle at step 1
v0 = np.zeros((m,3))   # j-th velocity of i-th particle at step 0
v1 = np.zeros((m,3))   # j-th velocity of i-th particle at step 1
Kch = N//nch
fvir2 = np.zeros(Kch)
up = np.zeros(Kch)
ub = np.zeros(Kch)
ut = np.zeros(Kch)
srdf = np.zeros((Ndf,2))

print("G =", G, ", m =", m, ", N =", N)
print("alpha =", alpha)
print()

@njit(parallel=True)
def forcal(x):
    R2 = R0*R0
    a = np.zeros((m,3))

    for i in prange(0,m):
        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        for j in prange(0,i):
            d1 = x[i,0] - x[j,0]
            d2 = x[i,1] - x[j,1]
            d3 = x[i,2] - x[j,2]
            r = sqrt(d1*d1 + d2*d2 + d3*d3)
            s3 = r*r*r
            c1 += d1/s3
            c2 += d2/s3
            c3 += d3/s3
        a[i,0] += c1
        a[i,1] += c2
        a[i,2] += c3

        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        for j in prange(i+1,m):
            d1 = x[i,0] - x[j,0]
            d2 = x[i,1] - x[j,1]
            d3 = x[i,2] - x[j,2]
            r = sqrt(d1*d1 + d2*d2 + d3*d3)
            s3 = r*r*r
            c1 += d1/s3
            c2 += d2/s3
            c3 += d3/s3
        a[i,0] += c1
        a[i,1] += c2
        a[i,2] += c3

    for i in prange(0,m):
        rk2 = x[i,0]*x[i,0] + x[i,1]*x[i,1] + x[i,2]*x[i,2]
        if (rk2 <= R2):
            a[i,0] -= x[i,0]
            a[i,1] -= x[i,1]
            a[i,2] -= x[i,2]
        else:
            rk32 = rk2**1.5
            a[i,0] -= m/rk32*x[i,0]
            a[i,1] -= m/rk32*x[i,1]
            a[i,2] -= m/rk32*x[i,2]

    return a

@njit(parallel=True)
def forin(x):
    R2 = R0*R0
    a = np.zeros((m,3))

    for i in prange(0,m):
        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        for j in prange(0,i):
            if (x[j,0]*x[j,0] + x[j,1]*x[j,1] + x[j,2]*x[j,2] < Rcore*Rcore) :
                d1 = x[i,0] - x[j,0]
                d2 = x[i,1] - x[j,1]
                d3 = x[i,2] - x[j,2]
                r = sqrt(d1*d1 + d2*d2 + d3*d3)
                s3 = r*r*r
                c1 += d1/s3
                c2 += d2/s3
                c3 += d3/s3
        a[i,0] += c1
        a[i,1] += c2
        a[i,2] += c3

        c1 = 0.0
        c2 = 0.0
        c3 = 0.0
        for j in prange(i+1,m):
            if (x[j,0]*x[j,0] + x[j,1]*x[j,1] + x[j,2]*x[j,2] < Rcore*Rcore) :
                d1 = x[i,0] - x[j,0]
                d2 = x[i,1] - x[j,1]
                d3 = x[i,2] - x[j,2]
                r = sqrt(d1*d1 + d2*d2 + d3*d3)
                s3 = r*r*r
                c1 += d1/s3
                c2 += d2/s3
                c3 += d3/s3
        a[i,0] += c1
        a[i,1] += c2
        a[i,2] += c3

    for i in prange(0,m):
        rk2 = x[i,0]*x[i,0] + x[i,1]*x[i,1] + x[i,2]*x[i,2]
        if (rk2 <= R2):
            a[i,0] -= x[i,0]
            a[i,1] -= x[i,1]
            a[i,2] -= x[i,2]
        else:
            rk32 = rk2**1.5
            a[i,0] -= m/rk32*x[i,0]
            a[i,1] -= m/rk32*x[i,1]
            a[i,2] -= m/rk32*x[i,2]

    return a

@njit(parallel=False)
def specular(x,v):
    w = v
    for k in range (0,m):
        dist2 = x[k,0]*x[k,0] + x[k,1]*x[k,1] + x[k,2]*x[k,2]
        vl = x[k,0]*w[k,0] + x[k,1]*w[k,1] + x[k,2]*w[k,2]
        v2 = w[k,0]*w[k,0] + w[k,1]*w[k,1] + w[k,2]*w[k,2]
        if (dist2 > Rb*Rb and vl > 0):
            w[k,0] -= 2/dist2*vl*x[k,0]
            w[k,1] -= 2/dist2*vl*x[k,1]
            w[k,2] -= 2/dist2*vl*x[k,2]
        if (v2 > 25/G):
            v0 = sqrt(v2)
            w[k,0] *= u0/v0
            w[k,1] *= u0/v0
            w[k,2] *= u0/v0
    return w

@njit(parallel=True)
def upb(x,v):
    R2 = R0*R0
    a = np.zeros(3)
    a0 = 0.0
    a1 = 0.0
    a2 = 0.0
    for i in prange(0,m-1):
        for j in prange(i+1,m):
            d1 = x[i,0] - x[j,0]
            d2 = x[i,1] - x[j,1]
            d3 = x[i,2] - x[j,2]
            r = sqrt(d1*d1 + d2*d2 + d3*d3)
            a0 += 1/r

    for k in prange(0,m):
        rk2 = x[k,0]*x[k,0] + x[k,1]*x[k,1] + x[k,2]*x[k,2]
        if (rk2 <= R2):
            a1 += rk2/2
        else:
            a1 += -m/rk2**0.5 + 1.5*R2

        a2 += (v[k,0]*v[k,0] + v[k,1]*v[k,1] + v[k,2]*v[k,2])/2

    a[0] = a0
    a[1] = a1
    a[2] = a2

    return a

@njit(parallel=False)
def rdfcalc(x):
    ardf = np.zeros((Ndf,2))
    rmin = 2*R0
    kmin = 0
    for k in range (0,m):
        dist = sqrt(x[k,0]*x[k,0] + x[k,1]*x[k,1] + x[k,2]*x[k,2])
        ki = int(dist/R0*Ndf)
        if (ki < Ndf):
            ardf[ki,0] = ardf[ki,0] + 1.0
        if (rmin > dist):
            rmin = dist
            kmin = k
    for k in range (0,m):
        if (k == kmin):
            continue
        dist = sqrt((x[k,0] - x[kmin,0])**2 + (x[k,1] - x[kmin,1])**2 + (x[k,2] - x[kmin,2])**2)
        ki = int(dist/R0*Ndf)
        if (ki < Ndf):
            ardf[ki,1] = ardf[ki,1] + 1.0

    return ardf

print(time.ctime())

# if swap exist
swap_file_name = f'swap_{m}_{idnum}.dat'
if os.path.exists(swap_file_name):
    # read swap
    la = np.loadtxt(swap_file_name, dtype=float)
    x001 = la[:,0]
    x002 = la[:,1]
    x003 = la[:,2]
    v001 = la[:,3]
    v002 = la[:,4]
    v003 = la[:,5]
    lx0 = len(x001)
    x0 = np.concatenate((np.reshape(x001, (lx0,1)), np.reshape(x002, (lx0,1)), np.reshape(x003, (lx0,1))), axis=1)
    v0 = np.concatenate((np.reshape(v001, (lx0,1)), np.reshape(v002, (lx0,1)), np.reshape(v003, (lx0,1))), axis=1)
else:
    # initialize
    i = 0
    while i < m:
        x0[i,0] = random.uniform(-Rb,Rb)
        x0[i,1] = random.uniform(-Rb,Rb)
        x0[i,2] = random.uniform(-Rb,Rb)
        if (x0[i,0])**2 + (x0[i,1])**2 + (x0[i,2])**2 > Rb*Rb:
            continue
        i += 1
    v0 = sqrt(1/G)*np.random.normal( 0.0 , 1.0, size = (m,3) )
    v0 = specular(x0,v0)

fc0 = forcal(x0)
fnd0 = delta*np.random.normal( 0.0 , 1.0, size = (m,3) )
f0 = -fr*v0 + fnd0 + fc0

print()
for i in range (1,N):
    x1 = x0 + v0*tau0 + tau0*tau0/2*f0
    fc1 = forcal(x1)
    fnd1 = delta*np.random.normal( 0.0 , 1.0, size = (m,3) )
    v1 = (v0 + tau0/2*(f0 + fnd1 + fc1))/(1 + gamma)
    f1 = -fr*v1 + fnd1 + fc1
    f0 = f1
    x0 = x1
    v0 = v1
    v0 = specular(x0,v0)
    if (i%nch == 0):
        # print(i//nch, " of ", Kch, " data collection")
        utmp = upb(x0,v0)
        up[i//nch - 1] = utmp[0]/m
        ub[i//nch - 1] = utmp[1]/m
        ut[i//nch - 1] = utmp[2]/m
        frn = forin(x0)
        fvir = 0.0
        mfv = 0
        for k in range (0,m):
            if (x0[k,0]*x0[k,0] + x0[k,1]*x0[k,1] + x0[k,2]*x0[k,2] < Rcore*Rcore):
                fvir += x0[k,0]*frn[k,0] + x0[k,1]*frn[k,1] + x0[k,2]*frn[k,2]
                mfv += 1
        fvir2[i//nch - 1] = fvir/mfv
        srdf = srdf + rdfcalc(x0)

fvir = 1/3*np.mean(fvir2)
pair = np.mean(up)
backgr = np.mean(ub)
temper = np.mean(ut)
pressure = (1/G + fvir)/(4*pi/3)

xcodf = 1/3**0.5 - 1
xrdf = (np.arange(1,Ndf+1) + xcodf)/Ndf*R0
srdf1 = srdf[:,0]/3/Kch/R0*Ndf/xrdf**2
srdf2 = srdf[:,1]*m/(m-1)/3/Kch/R0*Ndf/xrdf**2

for k in range (0,Ndf):
    print(xrdf[k],srdf1[k])

for k in range (0,Ndf):
    print(xrdf[k],srdf2[k])

print(time.ctime())

print(alpha, G*fvir, pressure)
print((pair - 2*backgr)*G/3, pair, backgr, 1.5/temper, pair + backgr - 0.9*m**(2/3))

# save swap
np.savetxt(swap_file_name, np.concatenate((np.reshape(x0,(len(x0),3)), np.reshape(v0,(len(x0),3))), axis=1))
