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
alpha = 0.5  
beta = 0.02  
delta = 0.01 
tau = 0.1    

x0 = 100 
y0 = 10 
t_total = 100  
n = 1000     

#Choose Method -> "lie_trotter", "strang", "third","sixth" ,"eight", "tenth", "twelfth" ,"fourteenth"
method="lie_trotter"
tol=1e-8 #tolerance for RK45 method


# -------------------------------
# Do not change below
# -------------------------------
  
dt = t_total / n
t_vals = np.linspace(0, t_total, n+1)

# -------------------------------
# Sigma and Gamma definitions for lokta-volterra
# -------------------------------
def sigma(y, t, x):
    return x * np.exp(t * (alpha - beta * y))

def gamma(x, t, y):
    return y * np.exp(t * (delta * x - tau))

x = np.zeros(n+1)
y = np.zeros(n+1)
x[0] = x0
y[0] = y0


#Complex coefficients for sixth method
def compute_a_k(k):
    # Compute the exponent for the complex number
    return np.exp((1j * np.pi) / (2 * k + 1)) /(2**(1 / (2 * k + 1)) + 2 * np.exp((1j * np.pi )/ (2 * k + 1)))

aa = compute_a_k(1) * compute_a_k(2)
bb = compute_a_k(2) * (1 - 2 * compute_a_k(1))
cc = compute_a_k(1) * (1 - 2 * compute_a_k(2))
dd = (1 - 2 * compute_a_k(1)) * (1 - 2 * compute_a_k(2))


#Complex coefficients for 8-14 method
def compute_a_k2(k):
    # Compute the exponent for the complex number
    a22= (1/4)+1j*(np.sin(np.pi / (2 * k + 1)))/(4 + 4 * np.cos(np.pi / (2 * k + 1)))
    return a22

a1=compute_a_k2(1)
a1b=np.conjugate(a1)
a2= compute_a_k2(2)
a2b=np.conjugate(a2)
a3=compute_a_k2(3)
a3b=np.conjugate(a3)
a4=compute_a_k2(4)
a4b=np.conjugate(a4)
a5=compute_a_k2(5)
a5b=np.conjugate(a5)
a6=compute_a_k2(6)
a6b=np.conjugate(a6)



if method == "lie_trotter":
    for k in range(n):
        x[k+1] = sigma(y[k], dt, x[k])
        y[k+1] = gamma(x[k+1], dt, y[k])

elif method == "strang":
    for k in range(n):
        x_k_half = sigma(y[k], dt/2, x[k])         
        y[k+1] = gamma(x_k_half, dt, y[k])         
        x[k+1] = sigma(y[k+1], dt/2, x_k_half)     
        
elif method == "third":
    aa=1/2 + 1j * np.sqrt(3)/6
    aa_bar=aa.conjugate()
    for k in range(n):
        y_1 = gamma(x[k], aa_bar * dt/2, y[k])
        x_2 = sigma(y_1, aa_bar * dt, x[k])
        y_3 = gamma(x_2, dt/2, y_1)
        x[k+1] = sigma(y_3, aa* dt, x_2)
        y[k+1] = gamma(x[k+1], aa * dt/2, y_3)

elif method == "sixth":
    for k in range(n):
            y_k_1 = gamma(x[k], aa * dt / 2, y[k])
            x_k_2 = sigma(y_k_1, aa * dt, x[k])
            y_k_3 = gamma(x_k_2, (aa + bb) * dt / 2, y_k_1)
            x_k_4 = sigma(y_k_3, bb * dt, x_k_2)
            y_k_5 = gamma(x_k_4, (aa + bb) * dt/2 , y_k_3)
            x_k_6 = sigma(y_k_5, aa * dt, x_k_4)
            y_k_7 = gamma(x_k_6, (aa + cc) * dt / 2, y_k_5)
            x_k_8 = sigma(y_k_7, cc * dt, x_k_6)
            y_k_9 = gamma(x_k_8, (cc + dd) * dt / 2, y_k_7)
            x_k_10 = sigma(y_k_9, dd * dt, x_k_8)
            y_k_11 = gamma(x_k_10, (cc + dd) * dt / 2, y_k_9)
            x_k_12 = sigma(y_k_11, cc * dt, x_k_10)
            y_k_13 = gamma(x_k_12, (aa + cc) * dt / 2, y_k_11)
            x_k_14 = sigma(y_k_13, aa * dt, x_k_12)
            y_k_15 = gamma(x_k_14, (aa + bb) * dt/2 , y_k_13)
            x_k_16 = sigma(y_k_15, bb * dt, x_k_14)
            y_k_17 = gamma(x_k_16, (aa + bb) * dt / 2, y_k_15)
            x[k+1] = sigma(y_k_17, aa * dt, x_k_16)
            y[k+1] = gamma(x[k+1], aa * dt / 2, y_k_17)
            

elif method == "eight":
    def z1loop(inpt, g,h):
            y_k_1 = gamma(g, inpt*a1 * dt/2,h)
            x_k_2 = sigma(y_k_1,inpt* a1* dt, g)
            y_k_3 = gamma(x_k_2,inpt*(1/4) * dt, y_k_1)
            x_k_4 = sigma(y_k_3, inpt*a1b * dt, x_k_2)
            y_k_5 = gamma(x_k_4, inpt*a1b * dt, y_k_3)
            x_k_6 = sigma(y_k_5, inpt*a1b * dt, x_k_4)
            y_k_7 = gamma(x_k_6,inpt*(1/4) * dt, y_k_5)
            xxx = sigma(y_k_7,inpt* a1 * dt, x_k_6)
            yyy = gamma(xxx, inpt*a1 *dt/2, y_k_7)
            return  xxx, yyy
    def zaloop6a(inpt, xp, yq):
            x1, y1=z1loop(inpt*a2, xp, yq)
            x2, y2=z1loop(inpt*a2b, x1, y1)
            x3, y3=z1loop(inpt*a2b, x2, y2)
            x4, y4=z1loop(inpt*a2,x3, y3)
            return x4, y4
    def zaloop8(xp, yq):
            x1, y1=zaloop6a(a3, xp, yq)
            x2, y2=zaloop6a(a3b, x1, y1)
            x3, y3=zaloop6a(a3b, x2, y2)
            x4, y4=zaloop6a(a3,x3, y3)
            return x4, y4
    for k in range(n):
            x[k+1], y[k+1]= zaloop8(x[k], y[k])
    

elif method == "tenth":
    def z1loop(inpt, g,h):
            y_k_1 = gamma(g, inpt*a1 * dt/2,h)
            x_k_2 = sigma(y_k_1,inpt* a1* dt, g)
            y_k_3 = gamma(x_k_2,inpt*(1/4) * dt, y_k_1)
            x_k_4 = sigma(y_k_3, inpt*a1b * dt, x_k_2)
            y_k_5 = gamma(x_k_4, inpt*a1b * dt, y_k_3)
            x_k_6 = sigma(y_k_5, inpt*a1b * dt, x_k_4)
            y_k_7 = gamma(x_k_6,inpt*(1/4) * dt, y_k_5)
            xxx = sigma(y_k_7,inpt* a1 * dt, x_k_6)
            yyy = gamma(xxx, inpt*a1 *dt/2, y_k_7)
            return  xxx, yyy
    def zaloop6a(inpt, xp, yq):
            x1, y1=z1loop(inpt*a2, xp, yq)
            x2, y2=z1loop(inpt*a2b, x1, y1)
            x3, y3=z1loop(inpt*a2b, x2, y2)
            x4, y4=z1loop(inpt*a2,x3, y3)
            return x4, y4
    def zaloop8a(inpt, xp, yq):
            x1, y1=zaloop6a(inpt*a3, xp, yq)
            x2, y2=zaloop6a(inpt*a3b, x1, y1)
            x3, y3=zaloop6a(inpt*a3b, x2, y2)
            x4, y4=zaloop6a(inpt*a3,x3, y3)
            return x4, y4
    def zaloop10( xp, yq):
            x1, y1=zaloop8a(a4, xp, yq)
            x2, y2=zaloop8a(a4b, x1, y1)
            x3, y3=zaloop8a(a4b, x2, y2)
            x4, y4=zaloop8a(a4,x3, y3)
            return x4, y4

    for k in range(n):
            x[k+1], y[k+1]= zaloop10(x[k], y[k])

elif method == "twelfth":
    def z1loop(inpt, g,h):
            y_k_1 = gamma(g, inpt*a1 * dt/2,h)
            x_k_2 = sigma(y_k_1,inpt* a1* dt, g)
            y_k_3 = gamma(x_k_2,inpt*(1/4) * dt, y_k_1)
            x_k_4 = sigma(y_k_3, inpt*a1b * dt, x_k_2)
            y_k_5 = gamma(x_k_4, inpt*a1b * dt, y_k_3)
            x_k_6 = sigma(y_k_5, inpt*a1b * dt, x_k_4)
            y_k_7 = gamma(x_k_6,inpt*(1/4) * dt, y_k_5)
            xxx = sigma(y_k_7,inpt* a1 * dt, x_k_6)
            yyy = gamma(xxx, inpt*a1 *dt/2, y_k_7)
            return  xxx, yyy
    def zaloop6a(inpt, xp, yq):
            x1, y1=z1loop(inpt*a2, xp, yq)
            x2, y2=z1loop(inpt*a2b, x1, y1)
            x3, y3=z1loop(inpt*a2b, x2, y2)
            x4, y4=z1loop(inpt*a2,x3, y3)
            return x4, y4
    def zaloop8a(inpt, xp, yq):
            x1, y1=zaloop6a(inpt*a3, xp, yq)
            x2, y2=zaloop6a(inpt*a3b, x1, y1)
            x3, y3=zaloop6a(inpt*a3b, x2, y2)
            x4, y4=zaloop6a(inpt*a3,x3, y3)
            return x4, y4
    def zaloop10a(inpt, xp, yq):
            x1, y1=zaloop8a(inpt*a4, xp, yq)
            x2, y2=zaloop8a(inpt*a4b, x1, y1)
            x3, y3=zaloop8a(inpt*a4b, x2, y2)
            x4, y4=zaloop8a(inpt*a4,x3, y3)
            return x4, y4
    def zaloop12(xp, yq):
            x1, y1=zaloop10a(a5, xp, yq)
            x2, y2=zaloop10a(a5b, x1, y1)
            x3, y3=zaloop10a(a5b, x2, y2)
            x4, y4=zaloop10a(a5,x3, y3)
            return x4, y4
    for k in range(n):
            x[k+1], y[k+1]= zaloop12(x[k], y[k])

elif method == "fourteenth":    
    def z1loop(inpt, g,h):
            y_k_1 = gamma(g, inpt*a1 * dt/2,h)
            x_k_2 = sigma(y_k_1,inpt* a1* dt, g)
            y_k_3 = gamma(x_k_2,inpt*(1/4) * dt, y_k_1)
            x_k_4 = sigma(y_k_3, inpt*a1b * dt, x_k_2)
            y_k_5 = gamma(x_k_4, inpt*a1b * dt, y_k_3)
            x_k_6 = sigma(y_k_5, inpt*a1b * dt, x_k_4)
            y_k_7 = gamma(x_k_6,inpt*(1/4) * dt, y_k_5)
            xxx = sigma(y_k_7,inpt* a1 * dt, x_k_6)
            yyy = gamma(xxx, inpt*a1 *dt/2, y_k_7)
            return  xxx, yyy
    def zaloop6a(inpt, xp, yq):
            x1, y1=z1loop(inpt*a2, xp, yq)
            x2, y2=z1loop(inpt*a2b, x1, y1)
            x3, y3=z1loop(inpt*a2b, x2, y2)
            x4, y4=z1loop(inpt*a2,x3, y3)
            return x4, y4
    def zaloop8a(inpt, xp, yq):
            x1, y1=zaloop6a(inpt*a3, xp, yq)
            x2, y2=zaloop6a(inpt*a3b, x1, y1)
            x3, y3=zaloop6a(inpt*a3b, x2, y2)
            x4, y4=zaloop6a(inpt*a3,x3, y3)
            return x4, y4
    def zaloop10a(inpt, xp, yq):
            x1, y1=zaloop8a(inpt*a4, xp, yq)
            x2, y2=zaloop8a(inpt*a4b, x1, y1)
            x3, y3=zaloop8a(inpt*a4b, x2, y2)
            x4, y4=zaloop8a(inpt*a4,x3, y3)
            return x4, y4
    def zaloop12a(inpt,xp, yq):
            x1, y1=zaloop10a(inpt*a5, xp, yq)
            x2, y2=zaloop10a(inpt*a5b, x1, y1)
            x3, y3=zaloop10a(inpt*a5b, x2, y2)
            x4, y4=zaloop10a(inpt*a5,x3, y3)
            return x4, y4
    def zaloop14(xp, yq):
            x1, y1=zaloop12a(a6, xp, yq)
            x2, y2=zaloop12a(a6b, x1, y1)
            x3, y3=zaloop12a(a6b, x2, y2)
            x4, y4=zaloop12a(a6,x3, y3)
            return x4, y4

    for k in range(n):
            x[k+1], y[k+1]= zaloop14(x[k], y[k])
      
    
else:
    raise ValueError("Invalid method. Choose one method.")

# -------------------------------
# RK45 Method for comparison
# -------------------------------
def system_of_equation(t, z):
    x, y = z
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - tau * y
    return [dxdt, dydt]

sol = solve_ivp(
    system_of_equation, 
    (0, t_total), 
    [x0, y0], 
    t_eval=t_vals, 
    rtol=tol, 
    atol=tol
)
x_rk45 = sol.y[0]
y_rk45 = sol.y[1]

# -------------------------------
# Plot x(t) and y(t)
# -------------------------------
plt.figure(figsize=(10, 6))
plt.plot(t_vals, x, 'r--', label=f'x(t) ({method})')
plt.plot(t_vals, x_rk45, 'r-', label='x(t) (RK45)')
plt.plot(t_vals, y, 'b--', label=f'y(t) ({method})')
plt.plot(t_vals, y_rk45, 'b-', label='y(t) (RK45)')
plt.xlabel('Time t')
plt.ylabel('Values of x and y')
plt.title(f'Lotka-Volterra: {method.title()} vs RK45')
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# Phase Plot: x vs y
# -------------------------------
plt.figure(figsize=(10, 6))
plt.plot(x, y, 'r--', label=f'{method.title()} Phase')
plt.plot(x_rk45, y_rk45, 'b-', label='RK45 Phase')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Phase Plot: {method.title()} vs RK45')
plt.legend()
plt.grid(True)
plt.show()

