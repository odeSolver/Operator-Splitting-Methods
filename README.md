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

## 📦 Requirements

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
| Model                        | Equation Form                                                                                                                                  | Parameters/Description                              |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| **Lotka-Volterra (2D)**      | $\frac{dx}{dt} = \alpha x - \beta x y$ <br> or <br> $\frac{dy}{dt} = \delta x y - \tau y$                                                | `alpha`, `beta`, `delta`, `tau` — predator-prey interaction |
| **Lorenz System**            | $\frac{dx}{dt} = \sigma(y - x)$<br> $\frac{dy}{dt} = x(\rho - z) - y$ <br> $\frac{dz}{dt} = xy - \beta z $                           | `sigma`, `rho`, `beta` — classic chaotic system             |
| **Van der Pol Oscillator**   | $\frac{d^2x}{dt^2} - \mu(1 - x^2)\frac{dx}{dt} + x = 0$ <br> OR <br> $\frac{dx}{dt} = y$ <br> $\frac{dy}{dt} = \mu(1 - x^2)y - x $ | `mu` — controls nonlinearity and damping                    |
| **Rossler Attractor**        | $\frac{dx}{dt} = -y - z$ <br> $\frac{dy}{dt} = x + \alpha y$ <br> $\frac{dz}{dt} = \beta - z(x - \gamma) $                                          | `alpha`, `beta`, `gamma` — models chaotic behavior                     |
| **Kermack-McKendrick (SIR)** | $\frac{dx}{dt} = -\alpha x y$ <br> $\frac{dy}{dt} = \alpha x y - \beta y$ <br> $\frac{z}{dt} = \beta y $                             | `alpha`, `beta` — disease transmission and recovery rates   |
| **Airy Equation**            | $\frac{d^2y}{dx^2} - x y = 0$                                                                                                                                 | No free parameters — fundamental ODE in quantum physics     |
| **3D Lokta-Volterra**            | $\frac{dx}{dt} = x(1-y)$ <br> $\frac{dy}{dt} = y(x-z)$ <br> $\frac{dz}{dt} =  z(y -1) $                                                                                                                                  | No free parameters    |
| ** Other**     | *Working on it*                                                                                                                           | NO parameters                     |

`````
``````

| **Model**                    | **Parameter** | **Description**                                     |
| ---------------------------- | ------------- | --------------------------------------------------- |
| **Lotka-Volterra**           | α             | Prey growth rate                                    |
|                              | β             | Predation rate coefficient                          |
|                              | δ             | Predator increase rate per consumed prey            |
|                              | τ             | Predator death rate                                 |
| **Lorenz System**            | σ             | Prandtl number (rate of thermal diffusion)          |
|                              | ρ             | Rayleigh number (related to temperature difference) |
|                              | β             | Geometric factor (affects vertical dissipation)     |
| **Van der Pol Oscillator**   | μ             | Nonlinearity and damping strength                   |
| **Rossler Attractor**        | a             | Linear growth in y                                  |
|                              | b             | Controls z feedback strength                        |
|                              | c             | Controls chaotic interaction between x and z        |
| **Kermack-McKendrick (SIR)** | β             | Infection rate (contact × transmission probability) |
|                              | γ             | Recovery or removal rate                            |
| **Airy Equation**            | –             | No parameters; uses $\frac{d^2y}{dx^2} - x y = 0$   |
| **Other / Custom Models**    | user-defined  | Input any differential equation and parameters      |


``````
``````

## 📊 Output

Each script may output:

  Solution curves as plots (saved as PNGs or JPGs in results/)

  Comparison graph with RK45 

  You can modify initial values to save or display these results as needed.

## 📌 Comparing RK45  

We implement the Runge-Kutta RK45 method for solving ODEs. You can find a detailed breakdown of the equations, step-by-step methodology, and an example problem in [Compare RK45](Compare_RK45.md).


## 📌 License

This project is licensed under the [MIT License](LICENSE). You’re free to use, modify, and distribute it with attribution.

## 🙋‍♂️ Author

**Arun**  
Ph.D. in Mathematics | M.S. in Applied Statistics  
📧 Reach out via the contact form on the website or email  
[LinkedIn](https://www.linkedin.com/in/arunbanjara/) 

## 🔗 Related Projects
If you’d like to see a web version of these solvers with interactivity and graph comparisons, visit the Flask Web App version
[ODE SOLVER](https://arun1111.pythonanywhere.com/)


