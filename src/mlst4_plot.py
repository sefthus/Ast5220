import sys
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt


#--------------------------- reading in data

x_hires, k_hires       = np.loadtxt('x_k_hires.dat',usecols=(0,1), unpack=True)
l_val       = np.loadtxt('l_val.dat',unpack=True)
#S   	    = np.loadtxt('source_func.dat', unpack=True)
Theta_l  = np.loadtxt('Theta_l.dat',unpack=True)
C_l_int  = np.loadtxt('C_l_integrand.dat', unpack=True)
C_l_int100  = np.loadtxt('C_l_integrand_l100.dat', unpack=True)

l_hires, C_l            = np.loadtxt('C_l.dat',usecols=(0,1), unpack=True)

planck_low  = np.loadtxt("COM_PowerSpect_CMB-TT-loL-full_R2.02.txt", unpack=True,skiprows=3)
planck_hi  = np.loadtxt("COM_PowerSpect_CMB-TT-hiL-full_R2.02.txt", unpack=True,skiprows=3)

planck_l    = np.hstack([planck_low[0], planck_hi[0]])
C_l_planck  = np.hstack([planck_low[1], planck_hi[1]])
error       = np.hstack([planck_low[2], planck_hi[2]])

#----------------------------- variables
Mpc     = 3.08568025e22 # Mpc to km
c       = 2.99792458e8   # m/s
h0      = 0.7
H0      = h0 * 100. * 1.e3 / Mpc
n_s=0.96
k_hires      = k_hires*c/H0
#sys.exit()
# Normalize C_l
#C_l = C_l/np.max(C_l)*5775
#------------------------------ plotting
#sys.exit()
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#-------------- plot integrands
'''
Sj  = np.loadtxt('Sj_l.dat', unpack=True)
bessel= np.loadtxt('besseltest.dat', unpack=True)
bessel2= np.loadtxt('besseltest2.dat', unpack=True)

plt.figure()

#plt.plot(x_hires, bessel, label='1')
#plt.plot(x_hires, bessel2,label='2')

plt.plot(x_hires, Sj/1e-3) 
#plt.legend(loc='upper left')
plt.ylabel(r'$\tilde S(k,x)j_l[\eta_0-\eta]/10^{-3}$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.show()
#sys.exit()

l_val=np.array([6, 100, 200, 500, 1000, 1200])
plt.figure()
[plt.plot(k_hires, l_val[i]*(1+l_val[i])*C_l_int[i]*H0/c, label=(r'$l=%d$'%l_val[i])) for i in range(len(l_val))]
#plt.plot(k_hires, l_val[2]*(1+l_val[2])*C_l_int[i]*H0/c)
plt.legend(loc='upper right')
plt.ylabel(r'$l(l+1)\Theta_l^2_l(k)H_0/ck$')
plt.xlabel(r'$kc/H_0$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()
#sys.exit()

# --------------- plot Theta_l and the power spectrum C_l
plt.figure()
[plt.plot(k_hires, Theta_l[i], label=(r'$l=%d$' %l_val[i])) for i in range(len(l_val))]
plt.legend(loc='upper right')
plt.ylabel(r'$\Theta_l$')
plt.xlabel(r'$kc/H_0$')
plt.tight_layout()
plt.show()

normal = 5775./np.max(C_l)
plt.figure()
plt.errorbar(planck_l, C_l_planck, yerr=error, label='Planck')
plt.plot(l_hires, C_l*normal, label='simulated')
plt.legend(loc='upper left')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlabel(r'$l$')
plt.tight_layout()
plt.show()
'''
#------------------------- parameter estimation ---------------------
C_l_h66           = np.loadtxt('C_l_h66.dat',usecols=(1), unpack=True)
C_l_h74            = np.loadtxt('C_l_h74.dat',usecols=(1), unpack=True)
C_l_n099           = np.loadtxt('C_l_ns099.dat',usecols=(1), unpack=True)
C_l_n092            = np.loadtxt('C_l_ns092.dat',usecols=(1), unpack=True)

C_l_b042            = np.loadtxt('C_l_om_b042.dat',usecols=(1), unpack=True)
C_l_b048            = np.loadtxt('C_l_om_b048.dat',usecols=(1), unpack=True)
C_l_m220            = np.loadtxt('C_l_om_m220.dat',usecols=(1), unpack=True)
C_l_m228            = np.loadtxt('C_l_om_m228.dat',usecols=(1), unpack=True)
C_l_r83d4            = np.loadtxt('C_l_om_r83d4.dat',usecols=(1), unpack=True)
C_l_r43d5            = np.loadtxt('C_l_om_r43d5.dat',usecols=(1), unpack=True)

C_l_m = np.array([C_l_m220,C_l_m228])
C_l_r = np.array([C_l_r83d4,C_l_r43d5])
C_l_b = np.array([C_l_b042,C_l_b048])
C_l_h = np.array([C_l_h66,C_l_h74])
C_l_n = np.array([C_l_n092,C_l_n099])

param = 0 
plt.figure()
plt.errorbar(planck_l, C_l_planck, yerr=error, label=r'Planck')
if param=0:
    plt.plot(l_hires, C_l_m[0]*5775./np.max(C_l_m[0]), label=(r'\Omega_m=0.220'))
    plt.plot(l_hires, C_l_m[1]*5775./np.max(C_l_m[1]), label=(r'\Omega_m=0.228'))
if param=1:
    plt.plot(l_hires, C_l_r[0]*5775./np.max(C_l_r[0]), label=(r'\Omega_r=4.3\cdot 10{-5}'))
    plt.plot(l_hires, C_l_r[1]*5775./np.max(C_l_r[1]), label=(r'\Omega_r=8.3\cdot 10{-4}'))
if param = 2:
    plt.plot(l_hires, C_l_b[0]*5775./np.max(C_l_b[0]), label=(r'\Omega_b=0.042'))
    plt.plot(l_hires, C_l_b[1]*5775./np.max(C_l_b[1]), label=(r'\Omega_b=0.048'))
if param = 3:
    plt.plot(l_hires, C_l_h[0]*5775./np.max(C_l_h[0]), label=(r'h=66'))
    plt.plot(l_hires, C_l_h[1]*5775./np.max(C_l_h[1]), label=(r'h=74'))
if param = 4:
    plt.plot(l_hires, C_l_n[0]*5775./np.max(C_l_n[0]), label=(r'n=0.92'))
    plt.plot(l_hires, C_l_n[1]*5775./np.max(C_l_n[1]), label=(r'h=0.99'))

plt.legend(loc='upper left')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlabel(r'$l$')
plt.tight_layout()
plt.show()

