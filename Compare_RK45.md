# 🧠 Runge-Kutta Method (RK45) for Solving Ordinary Differential Equations

The Runge-Kutta methods are a family of iterative techniques used to approximate solutions to ordinary differential equations (ODEs) of the form:

$$
\frac{dy}{dt} = f(t,y), \quad y(t_0) = y_0
$$

The RK45 method (also known as the Dormand–Prince method) is a popular adaptive step-size method that combines 4th and 5th order Runge-Kutta methods to estimate the local error and adjust the time step accordingly.

## 🌟 Main Idea
RK45 computes two estimates of the next value \( y_{n+1} \) using:

- A 4th-order method \(y_{n+1}^{[4]} \)

- A 5th-order method \( y_{n+1}^{[5]} \)

The difference between them provides an estimate of the local truncation error, which is used to adapt the step size \( h \).
## 🧮 RK45 Equations (Simplified Form)

At each step, RK45 computes intermediate stages:

$$
\begin{align*}
k_1 &= f(t_n, y_n) \\
k_2 &= f(t_n + c_2 h, \, y_n + a_{21} h k_1) \\
k_3 &= f(t_n + c_3 h, \, y_n + a_{31} h k_1 + a_{32} h k_2) \\
&\vdots \\
k_6 &= f\left(t_n + c_6 h, \, y_n + \sum_{j=1}^5 a_{6j} h k_j\right)
\end{align*}
$$

Then, the next values of \( y \) are approximated by:

$$
y_{n+1}^{[4]} = y_n + h \sum_{i=1}^6 b_i^{[4]} k_i
$$

$$
y_{n+1}^{[5]} = y_n + h \sum_{i=1}^6 b_i^{[5]} k_i
$$

The error estimate is:

$$
E = \left| y_{n+1}^{[5]} - y_{n+1}^{[4]} \right|
$$

---

## 📉 Step Size Control with Relative and Absolute Tolerances

The local error estimate \( E \) is compared against a user-defined tolerance \( \epsilon \) using the condition:

$$
\frac{|E|}{\text{ATOL} + \text{RTOL} \cdot |y_{n+1}^{[5]}|} < 1
$$

Where:

- **ATOL** = Absolute Tolerance  
- **RTOL** = Relative Tolerance  

This ensures the error is controlled relative to the size of the solution, which is crucial when values vary across magnitudes.

---

## ✅ If the Error is Acceptable

- Accept the step
- Increase the step size \( h \) to be more efficient

## ❌ If the Error is Too Large

- Reject the step
- Decrease \( h \) to improve accuracy

A common formula for adapting the step size is:

$$
h_{\text{new}} = h \cdot \min\left( \max\left( 0.1, \, 0.84 \cdot \left( \frac{\text{tolerance}}{E} \right)^{1/4} \right), 4.0 \right)
$$

---

## ✏️ Summary

- RK45 uses **two approximations** (4th and 5th order) to estimate and control local error.
- It is **adaptive**, meaning the step size changes to keep the error within a user-defined tolerance.
- Use **absolute tolerance** to limit error in small solutions.
- Use **relative tolerance** to scale error for large solutions.

## 🛠️ Example: Solving a 2D System with RK45

To demonstrate RK45, let's solve the following system of differential equations:

$$
\frac{dx}{dt} = f_1(x, y) = y
$$

$$
\frac{dy}{dt} = f_2(x, y) = -x
$$

This describes simple harmonic motion, similar to a pendulum or mass-spring system.

### Initial Conditions:
$$
x(0) = 1
$$
$$
y(0) = 0
$$
$$
h = 0.1
$$

### RK45 Step-by-Step Computation

We compute intermediate values:

**1. Calculate \( k_1 \):**

$$
k_{1x} = f_1(x_n, y_n) = y_n
$$
$$
k_{1y} = f_2(x_n, y_n) = -x_n
$$

**2. Calculate \( k_2 \):**

$$
k_{2x} = f_1(x_n + c_2 h, y_n + a_{21} h k_{1x})
$$
$$
k_{2y} = f_2(x_n + c_2 h, y_n + a_{21} h k_{1y})
$$

**3. Repeat for \( k_3, k_4, k_5, k_6 \)**  
Using similar expressions with appropriate coefficients.

**4. Compute Final Values:**

$$
x_{n+1}^{[4]} = x_n + h \sum_{i=1}^{6} b_i^{[4]} k_{ix}
$$
$$
y_{n+1}^{[4]} = y_n + h \sum_{i=1}^{6} b_i^{[4]} k_{iy}
$$

$$
x_{n+1}^{[5]} = x_n + h \sum_{i=1}^{6} b_i^{[5]} k_{ix}
$$
$$
y_{n+1}^{[5]} = y_n + h \sum_{i=1}^{6} b_i^{[5]} k_{iy}
$$

**5. Compute Error and Adjust Step Size:**

$$
E_x = |x_{n+1}^{[5]} - x_{n+1}^{[4]}|
$$
$$
E_y = |y_{n+1}^{[5]} - y_{n+1}^{[4]}|
$$

If:

$$
E_x < \text{tolerance}, \quad E_y < \text{tolerance}
$$

- Accept step, adjust \( h \)

Else:
- Reject step, decrease \( h \)

### ✨ Summary
This example shows how RK45 adapts the step size while solving a simple system of equations. In a real implementation, you'd loop through multiple steps until reaching the desired endpoint.
