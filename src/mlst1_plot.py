import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt


#------------------------------------- reading in data
omega_m, omega_b, omega_r, omega_lambda = np.loadtxt("omega_mbrl.dat",usecols=(0,1,2,3), unpack=True)
x_eta, eta  	= np.loadtxt("xeta_eta.dat",usecols=(0,1), unpack=True)
z, H 			= np.loadtxt("xeta_z_H.dat",usecols=(1,2), unpack=True)


#------------------------------------- defining stuffu
366
Mpc = 3.08568025e22
omega_total = omega_m+omega_b+omega_r+ omega_lambda
eta = eta/Mpc # m --> Mpc units
H = H* Mpc/1e3 # s^-1 --> km s^-1/Mpc


#------------------------------------- calculating radiation matter dark energy equality
indx1 = [i for i in range(len(omega_m)) if omega_r[i]>=omega_m[i]][-1]
x_radmat_eq = x_eta[indx1]
a_radmat_eq = np.exp(x_radmat_eq)
z_radmat_eq = 1./a_radmat_eq -1
omega_radmat_eq = omega_r[indx1]

indx2 = [i for i in range(len(omega_m)) if omega_m[i]>=omega_lambda[i]][-1]
x_matlam_eq = x_eta[indx2]
a_matlam_eq = np.exp(x_matlam_eq)
z_matlam_eq = 1./a_matlam_eq -1
omega_matlam_eq = omega_m[indx2]

print 'radiation - matter equality at:   x = %.3f , a = %.4f , z = %.3f, Omega_r = %.3f' % (x_radmat_eq, a_radmat_eq, z_radmat_eq, omega_radmat_eq)
print 'dark energy - matter equality at: x = %.3f , a = %.4f , z = %.3f, Omega_m = %.3f' % (x_matlam_eq, a_matlam_eq, a_matlam_eq, omega_matlam_eq)
print 'maximum dark matter density max(Omega_m) =', np.max(omega_m)
print 'final density parameter values: Omega_b = %.4f,  Omega_m = %.4f,  Omega_r = %.4f,  Omega_lambda = %.4f' \
			% (omega_b[-1], omega_m[-1], omega_r[-1], omega_lambda[-1])
print 'conformal time today eta(x=0) = ', eta[-1]/1e3, 'Gpc'

print 'Hubble parameter today H=', H[-1], 'km s^-1 Mpc^-1'
#------------------------------------- plotting
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':14})

#------------------------------------- plotting omega
plt.figure()
plt.plot( x_eta, omega_m, label = r'$\Omega_\mathrm{m}$')
plt.plot( x_eta, omega_b, label = r'$\Omega_\mathrm{b}$')
plt.plot( x_eta, omega_r, label = r'$\Omega_\mathrm{r}$')
plt.plot( x_eta, omega_lambda, label= r'$\Omega_\Lambda$')
plt.plot( x_eta, omega_total, label= r'$\Omega_\mathrm{Total}$')

plt.plot( x_radmat_eq, omega_radmat_eq, 'x', label = r'$\Omega_{r=m}$')
plt.plot( x_matlam_eq, omega_matlam_eq, 'x', label = r'$\Omega_{m=\Lambda}$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\Omega_\mathrm{X}$')
plt.legend(loc='center left')
plt.grid('on', linestyle ='--')
#plt.show()

#------------------------------------- plotting H against x
plt.figure()
plt.plot(x_eta,H)
plt.xlabel(r'$x$')
plt.ylabel(r'$H(x)$ [km s$^{-1}$Mpc$^{-1}$]')
plt.yscale('log')
plt.grid('on', linestyle ='--')

#plt.show()
#------------------------------------- plotting H against z
plt.figure()
plt.plot(z,H)
plt.gca().invert_xaxis() # because z
plt.xlabel(r'$z$')
plt.ylabel(r'$H(z)$ [km s$^{-1}$Mpc$^{-1}$]')
plt.yscale('log')
plt.xscale('log')
plt.grid('on', linestyle ='--')
plt.show()

#------------------------------------- plotting eta against x
plt.figure()
plt.plot(x_eta,eta)
plt.yscale('log')
plt.xlabel(r'$x$')
plt.ylabel(r'$\eta(x)$ [Mpc]')
plt.grid('on', linestyle ='--')
plt.tight_layout()
plt.show()



