"""
Created on Sat May 24 15:04:56 2025
@author: arunb
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# -------------------------------
# Initial Parameters (Choose your initial values and method)
# -------------------------------
sig=10
rho=28
beta=8/3 

x0 = 1.0  
y0 = 1.0  
z0 = 1.0
t_total = 20  
n = 2000     

#Choose Method -> "lie_trotter", "strang", "third", "sixth"
method="third"
tol=1e-8 #tolerance for RK45 method


# -------------------------------
# Do not change below
# -------------------------------
  
dt = t_total / n
t_vals = np.linspace(0, t_total, n+1)

# -------------------------------
# Sigma and Gamma definitions for lokta-volterra
# -------------------------------
def sigmaLor(t, x, y, z):
    return y * (1 - np.exp(-sig * t)) + x * np.exp(-sig * t)

def gammaLor(t, x, y, z):
    return x * (rho - z) * (1 - np.exp(-t)) + y * np.exp(-t)

def tauLor(t, x, y, z):
    return (x * y / beta) * (1 - np.exp(-beta * t)) + z * np.exp(-beta * t)

x = np.zeros(n+1)
y = np.zeros(n+1)
z= np.zeros(n + 1)
x[0] = x0
y[0] = y0
z[0]= z0

#complexcoefficient for sixth method
def compute_a_k(k):
    # Compute the exponent for the complex number
    return np.exp((1j * np.pi) / (2 * k + 1)) /(2**(1 / (2 * k + 1)) + 2 * np.exp((1j * np.pi )/ (2 * k + 1)))

aa = compute_a_k(1) * compute_a_k(2)
bb = compute_a_k(2) * (1 - 2 * compute_a_k(1))
cc = compute_a_k(1) * (1 - 2 * compute_a_k(2))
dd = (1 - 2 * compute_a_k(1)) * (1 - 2 * compute_a_k(2))

if method == "lie_trotter":
    for k in range(n):
        x[k+1] = sigmaLor(dt, x[k], y[k], z[k] )
        y[k+1] = gammaLor(dt, x[k+1], y[k], z[k])
        z[k+1] = tauLor(dt, x[k+1], y[k+1], z[k])

elif method == "strang":
    for k in range(n):
        z_1 = tauLor(dt/2, x[k], y[k], z[k])
        y_2 = gammaLor(dt/2, x[k], y[k], z_1)
        x[k+1] = sigmaLor(dt, x[k], y_2, z_1)
        y[k+1] = gammaLor(dt/2, x[k+1], y_2, z_1)
        z[k+1] = tauLor(dt/2, x[k+1], y[k+1], z_1)   
        
elif method == 'third':
        aT = 1/2 + 1j * np.sqrt(3)/6
        aTb=1/2 - 1j * np.sqrt(3)/6
        
        for k in range(n):
            y_1 = gammaLor(aT*dt/2, x[k], y[k], z[k])
            z_2 = tauLor(aT*aT*dt/2, x[k], y_1, z[k])
            x_3 = sigmaLor(aT*aT*dt, x[k], y_1, z_2)
            z_4 = tauLor(aT*dt/2, x_3, y_1, z_2)
            x_5 = sigmaLor(aT*aTb*dt, x_3, y_1, z_4 )
            z_6 = tauLor(aT*aTb*dt/2, x_5, y_1, z_4)
            y_7 = gammaLor(dt/2, x_5, y_1, z_6,)
            z_8 = tauLor(aTb*aT*dt/2, x_5, y_7, z_6,)
            x_9 = sigmaLor(aT*aTb*dt, x_5, y_7, z_8)
            z_10 = tauLor(aTb*dt/2, x_9, y_7, z_8)
            x[k+1] = sigmaLor(aTb*aTb*dt, x_9, y_7, z_10)
            z[k+1] = tauLor(aTb*aTb*dt/2, x[k+1], y_7, z_10)
            y[k+1] = gammaLor(aTb*dt/2, x[k+1], y_7, z[k+1])

elif method == 'sixth':
    def z1loop(inpt,a,b,c):
        z_1a = tauLor(inpt *aa * dt / 2, a, b, c)
        x_2a = sigmaLor(inpt *aa* dt, a, b, z_1a,)
        z_3a = tauLor(inpt  *(aa + bb) * dt / 2, x_2a,b, z_1a)
        x_4a = sigmaLor(inpt *bb * dt, x_2a, b, z_3a)
        z_5a = tauLor(inpt  *(aa + bb) * dt/2 , x_4a, b, z_3a)
        x_6a = sigmaLor(inpt  *aa * dt, x_4a,b, z_5a )
        z_7a = tauLor(inpt  *(aa + cc) * dt / 2, x_6a, b, z_5a)
        x_8a = sigmaLor(inpt  *cc * dt, x_6a, b, z_7a)
        z_9a = tauLor(inpt  *(cc + dd) * dt / 2, x_8a, b, z_7a)
        x_10a = sigmaLor(inpt  *dd * dt, x_8a,b, z_9a)
        z_11a = tauLor(inpt *(cc + dd) * dt / 2, x_10a,b, z_9a)
        x_12a = sigmaLor(inpt  *cc * dt, x_10a, b, z_11a)
        z_13a = tauLor(inpt  *(aa + cc) * dt / 2, x_12a,b, z_11a)
        x_14a = sigmaLor(inpt  *aa * dt, x_12a,b, z_13a )
        z_15a = tauLor(inpt  *(aa + bb) * dt/2, x_14a, b, z_13a)
        x_16a = sigmaLor(inpt  * bb * dt, x_14a, b, z_15a )
        z_17a = tauLor(inpt *(aa + bb) * dt / 2, x_16a, b, z_15a)
        xx = sigmaLor(inpt *aa * dt, x_16a, b, z_17a )
        zz = tauLor(inpt *aa * dt / 2, xx, b, z_17a)
        return xx , b, zz


    def z2loop(inpt,a,b,c):
        yy= gammaLor(inpt*dt/2, a, b, c)
        return a, yy, c

    def zaloop3( xp, yq, zq):
        x1, y1, z1=z2loop(aa, xp, yq, zq)
        x2, y2, z2=z1loop(aa, x1, y1, z1)
        x3, y3, z3=z2loop((aa + bb), x2, y2, z2)
        x4, y4, z4=z1loop(bb,x3, y3, z3)
        x5, y5, z5=z2loop((aa + bb), x4, y4, z4)
        x6, y6, z6=z1loop(aa,x5, y5, z5)
        x7, y7, z7=z2loop((aa +cc), x6, y6, z6)
        x8, y8, z8=z1loop(cc,x7, y7, z7)
        x9, y9, z9=z2loop((cc + dd), x8, y8, z8)
        x10, y10, z10=z1loop(dd,x9, y9, z9)
        x11, y11, z11=z2loop((cc + dd),x10, y10, z10)
        x12, y12, z12=z1loop(cc,x11, y11, z11)
        x13, y13, z13=z2loop((aa + cc),x12, y12, z12)
        x14, y14, z14=z1loop(aa, x13, y13, z13)
        x15, y15, z15=z2loop((aa + bb),x14, y14, z14)
        x16, y16, z16=z1loop(bb,x15, y15, z15)
        x17, y17, z17=z2loop((aa + bb),x16, y16, z16)
        x18, y18, z18=z1loop(aa, x17, y17, z17)
        x19, y19, z19=z2loop(aa,x18, y18, z18)

        return x19, y19, z19

    for k in range(n):
        x[k+1], y[k+1],z[k+1]= zaloop3(x[k], y[k],z[k])

       
    
else:
    raise ValueError("Invalid method. Choose one method.")

# -------------------------------
# RK45 Method for comparison
# -------------------------------
def system_of_equation(t, z):
    x, y, z = z
    dxdt = sig * (y - x)
    dydt = x * (rho  - z) - y
    dzdt = x*y - beta * z
    return [dxdt, dydt, dzdt]

sol = solve_ivp(system_of_equation,(0, t_total),[x0, y0, z0], t_eval=t_vals, rtol=tol, atol=tol)

x_rk45 = sol.y[0]
y_rk45 = sol.y[1]
z_rk45 = sol.y[2]

# -------------------------------
# Plot 
# -------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_rk45, y_rk45, z_rk45, lw=0.5, color='r', label='RK45')
ax.plot(x, y, z, lw=0.5, color='b', linestyle='--', label='Method')



ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Lorenz Attractor Comparison')
ax.legend()

plt.show()

# -------------------------------
# RMSE Calculation
# -------------------------------
rmse = np.sqrt(np.mean((x - x_rk45)**2 + (y - y_rk45)**2+ (z - z_rk45)**2))
print(f"\nRMSE between {method} and RK45: {rmse:.10f}")
