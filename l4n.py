import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

   # M_jupiter / (M_sun + M_jupiter)
G=6.67430e-11
M_sun =1.9884e30
M_jupiter = 1.89813e27
mu=M_jupiter/(M_sun+M_jupiter)
a=7.78412e11
n =np.sqrt(G*(M_sun+M_jupiter)/(a*a*a))#distance a=1 here-->a^3=1
x_sun = -mu*a
y_sun = 0.0
x_jupiter = (1 - mu)*a
y_jupiter = 0.0
x_L4 = (0.5 - mu)*a #(x_sun+x_jupiter)/2=x_L4)
y_L4 = np.sqrt(3) *a/ 2 #equalateral traingle

def equations(t, state):
    x, y, vx, vy = state
    # Distance from particle to Sun
    r_sun = np.sqrt( (x - x_sun)**2 + (y - y_sun)**2 )
    # Distance from particle to Jupiter
    r_jupiter = np.sqrt( (x - x_jupiter)**2 +(y - y_jupiter)**2)

    ax_sun = -G*M_sun * (x - x_sun) / r_sun**3#(x-x_sun)/r_sun x-component
    ay_sun = -G*M_sun * (y - y_sun) / r_sun**3# (y-y_sun)/r_sun y-component

    ax_jupiter = -G*M_jupiter * (x - x_jupiter) / r_jupiter**3#same here
    ay_jupiter = -G*M_jupiter * (y - y_jupiter) / r_jupiter**3

    ax_gravity = ax_sun + ax_jupiter
    ay_gravity = ay_sun + ay_jupiter

    ax = (
        ax_gravity+ 2 * n * vy+ n**2 * x
    )

    ay = (
        ay_gravity- 2 * n * vx+ n**2 * y
    )
    return [vx, vy, ax, ay]

delta = 1e6

x0 = x_L4+delta
y0 = y_L4

vx0 = 0.0
vy0 = 0.0
year=365.25*24*3600
initial_state = [x0,y0,vx0,vy0]
t_start = 0
t_end = 150*year

solution = solve_ivp(
    equations,
    [t_start, t_end],
    initial_state,
    rtol=1e-8,
    atol=1e-4,
    max_step=0.3*365.25*24*3600
)

x = solution.y[0]-x_L4##first element in state->x
y = solution.y[1]-y_L4##second element in state->y

plt.figure(figsize=(8, 8))

plt.plot(x, y, color='red',label="Test particle")

#plt.scatter(x_sun,y_sun,s=240,label="Sun")

#plt.scatter(x_jupiter,y_jupiter,s=80,label="Jupiter")

#plt.scatter(x_L4,y_L4,s=20,label="L4")
plt.xlabel("x")
plt.ylabel("y")

plt.title("L4 Libration")

plt.axis("equal")
plt.legend()
plt.grid()

plt.show()