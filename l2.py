import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

G = 6.67430e-11  
M_sun = 1.989e30  
M_earth = 5.972e24
a = 1.496e11  #distance from the earth to the sun
T = 365.25 * 24 * 3600 
omega = 2 * np.pi / T  ##angular speed of the Earth (also the a.s. of the spacecraft)
mu = M_earth / (M_sun + M_earth)##for calculating the barycentre of the Sun and Earth

def f(x):
    return -omega**2*x+G*M_sun/(x+mu*a)**2+G*M_earth/(x-(1-mu)*a)**2 ##the function

L2=newton(f,a+1e9)#find the root of the funciton: starting from x=a until f(x)=0

def g(x):
    return -omega**2*x+G*M_sun/(x+mu*a)**2-G*M_earth/(x-(1-mu)*a)**2

L1=newton(g,a+1e9)

##graphing##
t=np.linspace(0,T,1000)
x_sun=-mu*a*np.cos(omega*t)#the Sun is orbiting around the centre of mass with radius mu*a
y_sun=-mu*a*np.sin(omega*t)
x_earth=(1-mu)*a*np.cos(omega*t)#the Earth is orbiting around the barycentre 
y_earth=(1-mu)*a*np.sin(omega*t)
x_L2=L2*np.cos(omega*t)
y_L2=L2*np.sin(omega*t)
x_L1=L1*np.cos(omega*t)
y_L1=L1*np.sin(omega*t)

plt.figure(figsize=(10,10))
plt.plot(x_sun,y_sun,color='blue',label='Sun')
plt.plot(x_L2,y_L2,color='yellow',label='L2')
plt.plot(x_L1,y_L1,color='orange',label=L1)
plt.plot(x_earth, y_earth, label='Earth')
plt.scatter(0,0,color='black',label='Barycentre')#Barycentre is on the original point
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title('Sun-Earth-L2')
plt.legend()
plt.axis('equal')
plt.grid(True)
plt.show()