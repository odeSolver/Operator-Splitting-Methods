# Ordinary Differential Equation Solver (Python)

This repository contains standalone Python scripts that solve various systems of differential equations using a unified method based on my dissertation research.

## 📂 Examples

- **Lotka-Volterra** (`lotka_volterra.py`)  
  Predator-prey model with specified initial values.

- **Van der Pol Oscillator** (`van_der_pol.py`)  
  Nonlinear oscillator dynamics.

- **Lorenz System** (`lorenz.py`)  
  Classic chaotic system.

- **Airy Equation** (`example4_custom_model.py`)  
  A system using initial values of your choosing.

## 🧠 Methodology

Each example applies a technique that reduces N-dimensional differential systems into 1D separable parts, then solves them using exponential product formulas (Lie-Trotter, Strang, or Higher order splitting).  
This approach is from my graduate research.

## Requirements

Install the required packages using:

```bash
pip install -r requirements.txt
```

Contents of reeuirement.txt:
numpy, matplotlib, scipy

## ▶️ How to Run

```bash
python lotka_volterra.py

pyhom van_der_pol.py

python lorenz.py

python airy.py
```

## Output

Each script may output:

Solution curves as plots (saved as PNGs in results/)

Comparison graph with RK45 (if applicable)

You can modify each script or initial values to save or display these results as needed.


