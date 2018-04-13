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
b = np.where((x<-7.1)&(x>-7.5))[0]
c = np.argmin(abs(ddtau[b]))
idx2 = b[c]
idx3 = np.argmax(g)

print 'Last X_e value calculated by Saha: X_e=%.3f , z=%.2f, x=%.2f' %(Xe[idx], z[idx],x[idx])
print 'First X_e value calculated by Peebles: %.3f , z=%.2f' %(Xe[idx+1], z[idx+1])
print 'Free electron fraction today, X_e(z=0)=%.2e'  %Xe[-1]
print "turning point tau: x=%.2f, z=%.2f, tau=%.2f" %( x[idx2],z[idx2], tau[idx2])
print 'peak of visibility function g - recombination'
print '   x=%.2f, z=%.2f, g=%.2f, tau=%.2f, X_e=%.2f' % (x[idx3], z[idx3], g[idx3], tau[idx3],Xe[idx3])

#print Xe
#------------------------------ plotting
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#----------------------------------- plotting X_e against z
plt.figure()
plt.plot(z,Xe)

# some values
plt.plot(z[idx2],Xe[idx2],'x',label=r'$\tau$ decrease')
plt.plot(z[idx3],Xe[idx3],'x', label=r'$\tilde{g}_{max}$')

plt.gca().invert_xaxis() # because z
plt.xlabel(r'$z$')
plt.ylabel(r'$X_\mathrm{e}(z)$')
plt.yscale('log')
plt.xlim(2000,0)
plt.ylim(1e-4,1.4)
plt.grid('on', linestyle ='--')
plt.legend()
plt.tight_layout()
#plt.show()

#-------------------------------- plotting tau and derivatives
plt.figure()
plt.plot( x, tau, label = r'$\tau(x)$')
plt.plot( x, np.abs(dtau), '--', label = r"$|\tau'(x)|$")
plt.plot( x, np.abs(ddtau), '-.', label = r"$|\tau''(x)|$")
# some values
plt.plot(x[idx],tau[idx],'x',label='Peebles start')
plt.plot(x[idx2],tau[idx2],'x',label=r'$\tau$ decrease')
plt.plot(x[idx3],tau[idx3],'x', label=r'$\tilde{g}_{max}$')

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

# some values
plt.plot( x[idx2],g[idx2],'x',label=r'$\tau$ decrease')
plt.plot( x[idx3],g[idx3],'x', label=r'$\tilde{g}_{max}$')

plt.xlabel(r'$x$')
plt.ylabel(r"$ \tilde{g}(x),\,\tilde{g}(x)'/10,\,\tilde{g}(x)''/300$")
plt.legend(loc='upper right')
plt.xlim(-7.5,-6.2)
plt.grid('on', linestyle ='--')
plt.tight_layout()
plt.show()



