import numpy as np
import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt

#--------------------------- reading in data
k, ks       = np.loadtxt("k_ks.dat",usecols=(0,1), unpack=True)
x   	    = np.loadtxt("x_t.dat", unpack=True)
Phi_1,Phi_2,Phi_3,Phi_4,Phi_5,Phi_6	    = np.loadtxt("Phi.dat",usecols=(0,1,2,3,4,5), unpack=True)

Psi_1,Psi_2,Psi_3,Psi_4,Psi_5,Psi_6	    = np.loadtxt("Psi.dat",usecols=(0,1,2,3,4,5), unpack=True)

delta_1,delta_2,delta_3,delta_4,delta_5,delta_6	                = np.loadtxt("delta.dat",usecols=(0,1,2,3,4,5), unpack=True)
delta_b_1,delta_b_2,delta_b_3,delta_b_4,delta_b_5, delta_b_6    = np.loadtxt("delta_b.dat",usecols=(0,1,2,3,4,5), unpack=True)

v_1,v_2,v_3,v_4,v_5,v_6	             = np.loadtxt("v.dat",usecols=(0,1,2,3,4,5), unpack=True)
v_b_1,v_b_2,v_b_3,v_b_4,v_b_5,v_b_6	 = np.loadtxt("v_b.dat",usecols=(0,1,2,3,4,5), unpack=True)

Theta0_1,Theta0_2,Theta0_3,Theta0_4,Theta0_5,Theta0_6   = np.loadtxt("Theta0.dat",usecols=(0,1,2,3,4,5), unpack=True)
Theta1_1,Theta1_2,Theta1_3,Theta1_4,Theta1_5,Theta1_6   = np.loadtxt("Theta1.dat",usecols=(0,1,2,3,4,5), unpack=True)


#----------------------------- variables
Mpc     = 3.08568025e22 # Mpc to km
c       = 2.99792458e8   # m/s
h0      = 0.7
H0      = h0 * 100. * 1.e3 / Mpc
ks      = ks*c/H0

delta   = np.array([delta_1,delta_2,delta_3,delta_4,delta_5,delta_6])
delta_b = np.array([delta_b_1,delta_b_2,delta_b_3,delta_b_4,delta_b_5,delta_b_6])
v       = np.array([v_1,v_2,v_3,v_4,v_5,v_6])
v_b     = np.array([v_b_1,v_b_2,v_b_3,v_b_4,v_b_5,v_b_6])
Phi     = np.array([Phi_1,Phi_2,Phi_3,Phi_4,Phi_5,Phi_6])
Psi     = np.array([Psi_1,Psi_2,Psi_3,Psi_4,Psi_5,Psi_6])
Theta0  = np.array([Theta0_1,Theta0_2,Theta0_3,Theta0_4,Theta0_5,Theta0_6])
Theta1  = np.array([Theta1_1,Theta1_2,Theta1_3,Theta1_4,Theta1_5,Theta1_6])
#------------------------------ plotting

plt.rc('text',usetex=True)
plt.rc('font',**{'family':'serif','size':16})

#-------------- plot delta

plt.figure()
[plt.plot(x, delta[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\delta$')
plt.xlabel(r'$x$')
plt.yscale('log')
plt.tight_layout()
#plt.grid('on', linestyle ='--')

plt.figure()
[plt.plot(x, delta_b[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\delta_b$')
plt.xlabel(r'$x$')
plt.yscale('log')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()

# ------------- plot velocities
plt.figure()
[plt.plot(x, v[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$v$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
#plt.show()

plt.figure()
[plt.plot(x, v_b[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$v_b$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
plt.show()

# --------------- plot Phi and Psi
plt.figure()
[plt.plot(x, Phi[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\Phi$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')
#plt.show()

plt.figure()
[plt.plot(x, Psi[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\Psi$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')

# --------------- plot Theta_0 and Theta_1
plt.figure()
[plt.plot(x, Theta0[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\Theta_0$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')

plt.figure()
[plt.plot(x, Theta1[i], label=(r'$k=$'+"{0}".format(str(round(ks[i], 2))))) for i in range(len(k))]
plt.legend(loc='upper left')
plt.ylabel(r'$\Theta_1$')
plt.xlabel(r'$x$')
plt.tight_layout()
#plt.grid('on', linestyle ='--')

plt.show()

