import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('lowres_data.dat')
data_hr = np.loadtxt('highres_data.dat')

t = data[:,0]
x = data[:,1]
v = data[:,2]

t_hr = data_hr[:,0]
x_hr = data_hr[:,1]
v_hr = data_hr[:,2]


plt.title('x')
plt.plot(t, x,t_hr,x_hr)
plt.xlabel('t')
plt.ylabel('x')

plt.figure()
plt.title('v')
plt.plot(t, v,t_hr,v_hr)
plt.xlabel('t')
plt.ylabel('v')

plt.show()
'''
plt.figure()
plt.title('x_hr')
plt.plot(t_hr, x_hr)
plt.xlabel('t_hr')
plt.ylabel('x_hr')

plt.figure()
plt.title('v_hr')
plt.plot(t_hr, v_hr)
plt.xlabel('t_hr')
plt.ylabel('v_hr')

plt.show()'''
