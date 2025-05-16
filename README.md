# Differential Equation Solvers (Python)

This repository contains 4 standalone Python scripts that solve various systems of differential equations using a unified method based on my dissertation research.

## 📂 Examples

- **Lotka-Volterra** (`example1_lotka_volterra.py`)  
  Predator-prey model with specified initial values.

- **Van der Pol Oscillator** (`example2_van_der_pol.py`)  
  Nonlinear oscillator dynamics.

- **Lorenz System** (`example3_lorenz.py`)  
  Classic chaotic system.

- **Custom System** (`example4_custom_model.py`)  
  A system using initial values of your choosing.

## 🧠 Methodology

Each example applies a technique that reduces N-dimensional differential systems into 1D separable parts, then solves them using exponential product formulas (Lie-Trotter or Strang splitting).  
This approach is from my graduate research.

## ▶️ How to Run

```bash
python example1_lotka_volterra.py
