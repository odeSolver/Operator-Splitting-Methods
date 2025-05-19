# 📈 Compare Against RK45

This document explains how and why we compare the custom numerical method against the classical **Runge-Kutta 4(5) method (RK45)**.

---

## 🔧 Differential Equation Setup

We solve ODEs of the form:

\[
\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0
\]

RK45 estimates two values at each step:

- \( y_{n+1}^{[4]} \): 4th-order estimate  
- \( y_{n+1}^{[5]} \): 5th-order estimate

The difference gives a local error estimate:

\[
E = \left| y_{n+1}^{[5]} - y_{n+1}^{[4]} \right|
\]

---

## 🎯 Tolerance Control

RK45 uses two tolerances:

- `RTOL` (Relative Tolerance): handles large values
- `ATOL` (Absolute Tolerance): handles small values

A step is accepted if:

\[
\frac{|E|}{\text{ATOL} + \text{RTOL} \cdot |y_{n+1}^{[5]}|} < 1
\]

---

## 🔄 Adaptive Step Size

To maintain accuracy and efficiency, the step size is updated as:

\[
h_{\text{new}} = h \cdot \min\left( \max\left( 0.1, \, 0.84 \cdot \left( \frac{\text{tolerance}}{E} \right)^{1/4} \right), 4.0 \right)
\]

---

## ✅ Why Compare Against RK45?

- RK45 is a gold standard due to its **adaptive step size** and **error control**.
- It serves as a baseline to validate the **accuracy and stability** of custom solvers.
- Helps assess whether your fixed-step or alternative method performs within acceptable error bounds.

---

## 📊 Optional Error Metrics

To quantify comparison:

- **Relative Error**:
  \[
  \text{RelError} = \frac{|y_{\text{custom}} - y_{\text{RK45}}|}{|y_{\text{RK45}}| + \varepsilon}
  \]

- **Absolute Error**:
  \[
  \text{AbsError} = |y_{\text{custom}} - y_{\text{RK45}}|
  \]

Where \( \varepsilon \) is a small number to avoid division by zero.

---

### 📂 File generated as part of: `/Compare_RK45.md`
