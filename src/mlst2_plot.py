import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt

#--------------------------- reading in data
tau, dtau, ddtau		= np.loadtxt("x_tau.dat",usecols=(0,1,2), unpack=True)
g, dg, ddg  			= np.loadtxt("x_g.dat",usecols=(0,1,2), unpack=True)
x, z, Xe 				= np.loadtxt("x_z_Xe.dat",usecols=(0,1,2), unpack=True)

#-------------------------- print some values
idx = len(Xe[Xe>0.99])
#a = np.where((x>-10) &(x<-5))
b = np.where((x<-7.1)&(x>-7.5))[0]
c = np.argmin(abs(ddtau[b]))
d = np.argmax(abs(ddtau[b]))
idx2 = b[c]
idx3 =b[d]
print ddtau[idx2]
print x[idx2]
print 'Last X_e value calculated by Saha: %.3f , z=%.2f' %(Xe[idx], z[idx])
print 'First X_e value calculated by Peebles: %.3f , z=%.2f' %(Xe[idx+1], z[idx+1])
print 'Free electron fraction today, X_e(z=0)=%.2e'  %Xe[-1]

#print Xe
#------------------------------ plotting
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#----------------------------------- plotting X_e against z
plt.figure()
plt.plot(x,Xe)
#plt.gca().invert_xaxis() # because z
plt.xlabel(r'$z$')
plt.ylabel(r'$X_\mathrm{e}(z)$')
plt.yscale('log')
#plt.xlim(2000,0)
plt.ylim(1e-4,1.4)
plt.grid('on', linestyle ='--')
plt.tight_layout()
#plt.show()

#-------------------------------- plotting tau and derivatives
plt.figure()
plt.plot( x, tau, label = r'$\tau(x)$')
plt.plot( x, np.abs(dtau), '--', label = r"$|\tau'(x)|$")
plt.plot( x, np.abs(ddtau), '-.', label = r"$|\tau''(x)|$")
plt.plot(x[idx2],ddtau[idx2],'x')
plt.plot(x[idx2],tau[idx2],'x')
plt.plot(x[idx3],ddtau[idx3],'x')
plt.plot(x[idx3],tau[idx3],'x')
plt.xlabel(r'$x$')
plt.ylabel(r"$\tau(x),\,|\tau'(x)|,\,|\tau''(x)|$")
plt.legend(loc='upper right')
plt.yscale('log')
#plt.xlim(-19.5,0)
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
plt.show()



