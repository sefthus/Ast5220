import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt

#--------------------------- reading in data
tau, dtau, ddtau		= np.loadtxt("x_tau.dat",usecols=(0,1,2), unpack=True)
g, dg, ddg  			= np.loadtxt("x_g.dat",usecols=(0,1,2), unpack=True)
x, z, Xe 				= np.loadtxt("x_z_Xe.dat",usecols=(0,1,2), unpack=True)

#------------------------------ plotting
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#-------------------------------- plotting tau and derivatives
plt.figure()
plt.plot( x, tau, label = r'$\tau(x)$')
plt.plot( x, np.abs(dtau), '--', label = r"$|\tau'|$")
plt.plot( x, np.abs(ddtau), '-.', label = r"$|\tau''|$")
plt.xlabel(r'$x$')
plt.ylabel(r"$\tau(x),\,|\tau'(x)|,\,|\tau''(x)|$")
plt.legend(loc='upper right')
plt.yscale('log')
#plt.xlim(-10,0)
plt.grid('on', linestyle ='--')
plt.tight_layout()
#plt.show()


#--------------------------------- plotting g and derivatives
plt.figure()
plt.plot( x, g, label = r'$\tilde{g}(x)$')
plt.plot( x, dg/10., '--', label = r"$\tilde{g}'(x)$")
plt.plot( x, ddg/300., '-.', label = r"$\tilde{g}''(x)$")
plt.xlabel(r'$x$')
plt.ylabel(r"$ \tilde{g}(x),\,\tilde{g}(x)'/10,\,\tilde{g}(x)''/300$")
plt.legend(loc='upper right')
plt.xlim(-7.5,-6.2)
plt.grid('on', linestyle ='--')
plt.tight_layout()
#plt.show()

#----------------------------------- plotting X_e against z
plt.figure()
plt.plot(z,Xe)
plt.gca().invert_xaxis() # because z
plt.xlabel(r'$z$')
plt.ylabel(r'$X_\mathrm{e}(z)$')
plt.yscale('log')
plt.xlim(2000,0)
plt.ylim(1e-4,1.4)
plt.grid('on', linestyle ='--')
plt.tight_layout()
plt.show()


