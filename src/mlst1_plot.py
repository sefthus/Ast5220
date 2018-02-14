import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt

Mpc = 3.08568025e22

omega_m, omega_b, omega_r, omega_lambda = np.loadtxt("omega_mbrl.dat",usecols=(0,1,2,3), unpack=True)
x_eta, eta = np.loadtxt("xeta_eta.dat",usecols=(0,1), unpack=True)
x_t, eta_t = np.loadtxt("xt_eta_t.dat",usecols=(0,1), unpack=True)
z, H = np.loadtxt("xeta_z_H.dat",usecols=(1,2), unpack=True)

eta = eta/Mpc # m --> Mpc units
eta_t = eta_t/Mpc
H = H* Mpc/1e2 # s^-1 --> km s^-1/Mpc

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':14})

plt.figure()
plt.plot( x_eta, omega_m, label = r'$\Omega_\mathrm{m}$')
plt.plot( x_eta, omega_b, label = r'$\Omega_\mathrm{b}$')
plt.plot( x_eta, omega_r, label = r'$\Omega_\mathrm{r}$')
plt.plot( x_eta, omega_lambda, label= '$\Omega_\Lambda$')
plt.plot( x_eta, omega_m+omega_b+omega_r+ omega_lambda, label= '$\Omega_\mathrm{Total}$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\Omega_\mathrm{X}$')
plt.legend(loc='center left')
plt.grid('on', linestyle ='--')
#plt.show()

plt.figure()
plt.plot(x_eta,H)
plt.xlabel(r'$x$')
plt.ylabel(r'$H(x)$ [km s$^{-1}$Mpc$^{-1}$]')
plt.yscale('log')
plt.grid('on', linestyle ='--')

#plt.show()

plt.figure()
plt.plot(z,H)
plt.gca().invert_xaxis() # because z
plt.xlabel(r'$z$')
plt.ylabel(r'$H(z)$ [km s$^{-1}$Mpc$^{-1}$]')
plt.yscale('log')
plt.xscale('log')
plt.grid('on', linestyle ='--')
plt.show()

plt.figure()
plt.plot(x_eta,eta, label=r'$\eta_{computed}$')
plt.plot(x_t,eta_t,'--', label=r'$\eta_{splined}$')
plt.legend(loc='lower right')
plt.yscale('log')
plt.xlabel(r'$x$')
plt.ylabel(r'$\eta(x)$ [Mpc]')
plt.grid('on', linestyle ='--')
plt.show()


