

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

x0 = -.2 # Initial value for x
y0 = 0 # Initial value for y

# Define time and number of iterations
t_total =100  # Total time
n = 1000 # Number of iterations


#Choose Method -> "Lie-Trotter", "Strang", "Third","6th" ,"8th", "10th", "12th" ,"14th"
method="14th"
tol=1e-16 #tolerance for RK45 method



# Time step
dt =t_total / n
t_vals = np.linspace(0, t_total, n+1)


def aat(x):
    return 1-(x**2)
def bat(x):
    return -x


def sigmaV(y, t, x):
    return y * t +x

def gammaV(x, t, y):
    return bat(x)*( np.exp(t *aat(x))-1)/(aat(x)) + y* np.exp(t *aat(x))# Define initial values for x and y

# Initialize arrays to store the solutions
x = np.zeros(n+1)
y = np.zeros(n+1)


# Set initial values
x[0] = x0
y[0] = y0




def compute_a_k(k):
    # Compute the exponent for the complex number
    return np.exp((1j * np.pi) / (2 * k + 1)) /(2**(1 / (2 * k + 1)) + 2 * np.exp((1j * np.pi )/ (2 * k + 1)))

aa=compute_a_k(1) * compute_a_k(2)
bb=compute_a_k(2)*(1-2*compute_a_k(1))
cc=compute_a_k(1)*(1-2*compute_a_k(2))
dd=(1-2*compute_a_k(1))*(1-2*compute_a_k(2))

#Sixth to 14th definition for coefficients
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




if method == 'Lie-Trotter':
    for k in range(n):
        x[k+1] = sigmaV(y[k], dt, x[k])
        y[k+1] = gammaV(x[k+1], dt, y[k])

elif method == 'Strang':

    for k in range(n):
        x_k_half = sigmaV(y[k], dt / 2, x[k])
        y[k+1] = gammaV(x_k_half, dt, y[k])
        x[k+1] = sigmaV(y[k+1], dt / 2, x_k_half)
elif method == 'Third':
    aT = 1/2 + 1j * np.sqrt(3)/6
    aTb=1/2 - 1j * np.sqrt(3)/6
    for k in range(n):
        y_1 = gammaV(x[k], aT*dt/2, y[k])
        x_2 = sigmaV(y_1, aT* dt, x[k])
        y_3 = gammaV(x_2, dt/2, y_1)
        x[k+1] = sigmaV(y_3, aTb*dt, x_2)
        y[k+1] = gammaV(x[k+1], aTb*dt/2, y_3)

elif method == '8th':
    def z1loop(inpt, g,h):
        y_k_1 = gammaV(g, inpt*a1 * dt/2,h)
        x_k_2 = sigmaV(y_k_1,inpt* a1* dt, g)
        y_k_3 = gammaV(x_k_2,inpt*(1/4) * dt, y_k_1)
        x_k_4 = sigmaV(y_k_3, inpt*a1b * dt, x_k_2)
        y_k_5 = gammaV(x_k_4, inpt*a1b * dt, y_k_3)
        x_k_6 = sigmaV(y_k_5, inpt*a1b * dt, x_k_4)
        y_k_7 = gammaV(x_k_6,inpt*(1/4) * dt, y_k_5)
        xxx = sigmaV(y_k_7,inpt* a1 * dt, x_k_6)
        yyy = gammaV(xxx, inpt*a1 *dt/2, y_k_7)
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

elif method == '10th':
    def z1loop(inpt, g,h):
        y_k_1 = gammaV(g, inpt*a1 * dt/2,h)
        x_k_2 = sigmaV(y_k_1,inpt* a1* dt, g)
        y_k_3 = gammaV(x_k_2,inpt*(1/4) * dt, y_k_1)
        x_k_4 = sigmaV(y_k_3, inpt*a1b * dt, x_k_2)
        y_k_5 = gammaV(x_k_4, inpt*a1b * dt, y_k_3)
        x_k_6 = sigmaV(y_k_5, inpt*a1b * dt, x_k_4)
        y_k_7 = gammaV(x_k_6,inpt*(1/4) * dt, y_k_5)
        xxx = sigmaV(y_k_7,inpt* a1 * dt, x_k_6)
        yyy = gammaV(xxx, inpt*a1 *dt/2, y_k_7)
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

elif method == '12th':
    def z1loop(inpt, g,h):
        y_k_1 = gammaV(g, inpt*a1 * dt/2,h)
        x_k_2 = sigmaV(y_k_1,inpt* a1* dt, g)
        y_k_3 = gammaV(x_k_2,inpt*(1/4) * dt, y_k_1)
        x_k_4 = sigmaV(y_k_3, inpt*a1b * dt, x_k_2)
        y_k_5 = gammaV(x_k_4, inpt*a1b * dt, y_k_3)
        x_k_6 = sigmaV(y_k_5, inpt*a1b * dt, x_k_4)
        y_k_7 = gammaV(x_k_6,inpt*(1/4) * dt, y_k_5)
        xxx = sigmaV(y_k_7,inpt* a1 * dt, x_k_6)
        yyy = gammaV(xxx, inpt*a1 *dt/2, y_k_7)
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


elif method == '14th':
    def z1loop(inpt, g,h):
        y_k_1 = gammaV(g, inpt*a1 * dt/2,h)
        x_k_2 = sigmaV(y_k_1,inpt* a1* dt, g)
        y_k_3 = gammaV(x_k_2,inpt*(1/4) * dt, y_k_1)
        x_k_4 = sigmaV(y_k_3, inpt*a1b * dt, x_k_2)
        y_k_5 = gammaV(x_k_4, inpt*a1b * dt, y_k_3)
        x_k_6 = sigmaV(y_k_5, inpt*a1b * dt, x_k_4)
        y_k_7 = gammaV(x_k_6,inpt*(1/4) * dt, y_k_5)
        xxx = sigmaV(y_k_7,inpt* a1 * dt, x_k_6)
        yyy = gammaV(xxx, inpt*a1 *dt/2, y_k_7)
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


elif method == '6th':
    for k in range(n):
        y_k_1 = gammaV(x[k], aa * dt / 2, y[k])
        x_k_2 = sigmaV(y_k_1, aa * dt, x[k])
        y_k_3 = gammaV(x_k_2, (aa + bb) * dt / 2, y_k_1)
        x_k_4 = sigmaV(y_k_3, bb * dt, x_k_2)
        y_k_5 = gammaV(x_k_4, (aa + bb) * dt/2 , y_k_3)
        x_k_6 = sigmaV(y_k_5, aa * dt, x_k_4)
        y_k_7 = gammaV(x_k_6, (aa + cc) * dt / 2, y_k_5)
        x_k_8 = sigmaV(y_k_7, cc * dt, x_k_6)
        y_k_9 = gammaV(x_k_8, (cc + dd) * dt / 2, y_k_7)
        x_k_10 = sigmaV(y_k_9, dd * dt, x_k_8)
        y_k_11 = gammaV(x_k_10, (cc + dd) * dt / 2, y_k_9)
        x_k_12 = sigmaV(y_k_11, cc * dt, x_k_10)
        y_k_13 = gammaV(x_k_12, (aa + cc) * dt / 2, y_k_11)
        x_k_14 = sigmaV(y_k_13, aa * dt, x_k_12)
        y_k_15 = gammaV(x_k_14, (aa + bb) * dt/2 , y_k_13)
        x_k_16 = sigmaV(y_k_15, bb * dt, x_k_14)
        y_k_17 = gammaV(x_k_16, (aa + bb) * dt / 2, y_k_15)
        x[k+1] = sigmaV(y_k_17, aa * dt, x_k_16)
        y[k+1] = gammaV(x[k+1], aa * dt / 2, y_k_17)


else:
    raise ValueError("Invalid method. Choose one method.")

# Differential equations
def system_of_eqs(t, z):
    x, y = z
    dxdt = y
    dydt = aat(x)*y + bat(x)
    return [dxdt, dydt]


# Solve the system using solve_ivp
t_span = (0, t_total)
t_eval = np.linspace(0, t_total, n+1)
sol = solve_ivp(system_of_eqs, t_span, [x0, y0], method="RK45", t_eval=t_eval, rtol=tol, atol=tol)

# Extract the solutions
x_ode = sol.y[0]
y_ode = sol.y[1]





plt.figure(figsize=(10, 6))

# Phase plot for iterative method
plt.plot(x, y, 'r--', label='Method Phase Plot')

#Phase plot for ODE RK45 solution
plt.plot(x_ode, y_ode, 'b-', label='RK45 ODE Solution Phase Plot')

plt.title('Phase Plot: x vs. y')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()


rmse = np.sqrt(np.mean((x - x_ode)**2 + (y - y_ode)**2))
print(f"\nRMSE between {method} and RK45: {rmse:.20f}")
