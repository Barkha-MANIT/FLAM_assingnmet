# Parametric Curve Fitting — Parameter Estimation Assignment

This repository contains the solution for estimating the parameters  
θ, M, and X of the given parametric curve using 1500 observed (x, y) points.

---

## Model Definition

### **x-equation**

x(t) = t*cos(θ) - exp(M*|t|)*sin(0.3*t)*sin(θ) + X  

### **y-equation**

y(t) = 42 + t*sin(θ) + exp(M*|t|)*sin(0.3*t)*cos(θ)

Parameter bounds:

- \(0^\circ < \theta < 50^\circ\)  
- \(-0.05 < M < 0.05\)  
- \(0 < X < 100\)  
- \(6 \le t \le 60\)

---

# Final Estimated Parameters

| Parameter | Value |
|----------|--------|
| **θ (deg)** | **29.99997248** |
| **θ (rad)** | **0.52359830** |
| **M** | **0.02999999** |
| **X** | **54.99999819** |

---

# Final Parametric Equation

### x(t):
x(t) = t*cos(0.52359830) − exp(0.02999999*abs(t)) * sin(0.3*t) * sin(0.52359830) + 54.99999819 

### y(t):
y(t) = 42 + t*sin(0.52359830) + exp(0.02999999*abs(t)) * sin(0.3*t) * cos(0.52359830)

---

# Residual Summary

| Metric | Value |
|--------|--------|
| **Sum L1 Residual** | 0.003266 |
| **Mean Residual**   | 0.000002 |
| **Max Residual**    | 0.000011 |

These extremely small values indicate an almost perfect match.

---

# Method Overview

### 1. Initial Linear Approximation
A simplified model (ignoring the oscillatory exponential term) provides initial estimates for θ and X.

### 2. Global Search (Differential Evolution)
A subsample of 200 points is used for efficient global search.

### 3. Local Refinement (L-BFGS-B)
The DE output is refined using the full dataset of 1500 points.

### 4. Final Projection
Each observed point is projected onto the final fitted curve by optimizing over \(t \in [6, 60]\).

---

## How to Run the Code

Install dependencies:

```
pip install numpy scipy pandas
```

Run:

```
python fit_curve.py
```

---

## Repository Structure

| File | Description |
|------|-------------|
| fit_curve.py      | Main parameter estimation script |
| xy_data.csv       | Input dataset |
| fit_results.csv   | Projected t-values and residuals |
| README.md         | Assignment documentation |

---

## Final Answer (Plain Text Equation)

(t*cos(0.52359830) − exp(0.02999999*abs(t))*sin(0.3*t)*sin(0.52359830) + 54.99999819, 42 + t*sin(0.52359830) + exp(0.02999999*abs(t))*sin(0.3*t)*cos(0.52359830))

## Final Answer (LaTex Equation)

\left(t*\cos(0.52359830)-e^{0.02999999\left|t\right|}\cdot\sin(0.3t)\sin(0.52359830)\ +54.99999819,\;42+\ t*\sin(0.52359830)+e^{0.02999999\left|t\right|}\cdot\sin(0.3t)\cos(0.52359830)\right)

## Final Output Plot

![Fitted Curve Plot](fitted_curve_plot.png)

## Contact
Barkha Verma, devbarkha07@gmail.com

