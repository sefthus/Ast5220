import sys
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt


#--------------------------- reading in data
x_hires, k_hires       = np.loadtxt('x_k_hires.dat',usecols=(0,1), unpack=True)
k_val       = np.loadtxt('k_val_cl.dat',unpack=True)
#S   	    = np.loadtxt('source_func.dat', unpack=True)
Theta_l  = np.loadtxt('Theta_l.dat',unpack=True)

Theta_l_int = np.loadtxt('Theta_l_integrand.dat', unpack=True)
C_l_int  = np.loadtxt('C_l_integrand.dat', unpack=True)
C_l_int100 = np.loadtxt('C_l_integrand_l100.dat', unpack=True)

ls	    = np.loadtxt('ls.dat', unpack=True)
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

k_hires      = k_hires*c/H0
k_val      = k_val*c/H0
#sys.exit()
# Normalize C_l
#C_l = C_l/np.max(C_l)*5775
#------------------------------ plotting
#sys.exit()
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#-------------- plot integrands

Sj  = np.loadtxt('Sj_l.dat', unpack=True)
bessel= np.loadtxt('besseltest.dat', unpack=True)
bessel2= np.loadtxt('besseltest2.dat', unpack=True)
dPhi= np.loadtxt('dPhi.dat', unpack=True)
dPsi= np.loadtxt('dPsi.dat', unpack=True)
x_t= np.loadtxt('x_t.dat', unpack=True)

[plt.plot(x_t, dPhi[k], label=k) for k in range(6)] 
plt.legend(loc='upper left')
plt.ylabel(r'd$\Psi$')
plt.xlabel(r'$x_t$')
plt.tight_layout()
plt.show()
sys.exit()

'''
plt.figure()
#plt.plot(x_hires, Theta_l_int/1e-3) 

#plt.plot(x_hires, bessel, label='1')
#plt.plot(x_hires, bessel2,label='2')
plt.plot(x_hires, Sj/1e-3) 
plt.legend(loc='upper left')
plt.ylabel(r'$\tilde S(k,x)j_l[\eta_0-\eta]/10^{-3}$')
plt.xlabel(r'$x$')
plt.tight_layout()
plt.show()
#sys.exit()

n_s = 0.96
# l=100
plt.figure()
plt.plot(k_hires, C_l_int100/(c*k_hires/H0)**(n_s-1.))
#plt.legend(loc='upper left')
plt.ylabel(r'$\Theta_l^2_l(k)H_0/ck$')
plt.xlabel(r'$kc/H_0$')
plt.tight_layout()
plt.show()
#sys.exit()

plt.figure()
[plt.plot(l_hires, C_l_int[i]/(c*k_val[i]/H0)**(n_s-1.), label=(r'$kc/H_0=$'+"{0}".format(str(round(k_val[i], 2))))) for i in range(len(k_val))]
#plt.legend(loc='upper left')
plt.ylabel(r'$\Theta_l^2_l(k)H_0/ck$')
plt.xlabel(r'$l$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()

# --------------- plot Theta_l and the power spectrum C_l
plt.figure()
[plt.plot(l_hires, Theta_l[i], label=(r'$kc/H_0=$'+"{0}".format(str(round(k_val[i], 2))))) for i in range(len(k_val))]
plt.legend(loc='lower left')
plt.ylabel(r'$\Theta_l$')
plt.xlabel(r'$l$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()
'''
normal = 5775/np.max(C_l)
plt.figure()

plt.errorbar(planck_l, C_l_planck, yerr=error, label='Planck')
plt.plot(l_hires, C_l*normal, label='simulated')
plt.legend(loc='upper left')
plt.ylabel(r'$l(l+1)C_l/2\pi$')
plt.xlabel(r'$l$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()

