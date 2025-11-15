#!/usr/bin/env python3
"""
fit_curve_submission.py

Find parameters theta (degrees), M, X for the parametric model:
    x(t) = t*cos(theta) - exp(M*|t|)*sin(0.3*t)*sin(theta) + X
    y(t) = 42 + t*sin(theta) + exp(M*|t|)*sin(0.3*t)*cos(theta)

Input:
    - xy_data.csv  (must have columns 'x' and 'y')

Output:
    - fit_results.csv  (columns: x,y,t_proj,residual)
    - README_submission.txt
    - printed parameter estimates and residual statistics

Notes:
    - The script first obtains a linear heuristic to initialize parameters,
      then runs a global search on a subsample and a local refinement on full data.
    - Designed to be robust and reasonably fast on typical laptops.
"""

import os
import time
import numpy as np
import pandas as pd
from math import atan, degrees, radians
from scipy.optimize import differential_evolution, minimize, minimize_scalar
from functools import partial

# -------------------------
# User-tunable configuration
# -------------------------
DATA_FILENAME = "xy_data.csv"
OUT_CSV = "fit_results.csv"
OUT_README = "README_submission.txt"

SUBSAMPLE_SIZE = 200       # number of points used for global search (DE)
DE_MAXITER = 60            # DE iterations for subsample stage
DE_POPSZ = 10              # DE population size
DE_WORKERS = 1             # use single worker for broad compatibility (Windows)
INNER_MAXITER = 180        # inner scalar minimizer iterations during optimization
FINAL_INNER_MAXITER = 400  # final inner iterations for projection
INNER_XATOL = 1e-4
FINAL_XATOL = 1e-6
T_MIN, T_MAX = 6.0, 60.0

# variable bounds (theta in degrees)
BOUNDS = [(1e-6, 50.0), (-0.05, 0.05), (0.0, 100.0)]

RANDOM_SEED = 42

# -------------------------
# Model definition
# -------------------------
def model_xy(t, theta_rad, M, X):
    """Return (x, y) for input scalar or array t."""
    t = np.array(t, dtype=float)
    x = t * np.cos(theta_rad) - np.exp(M * np.abs(t)) * np.sin(0.3 * t) * np.sin(theta_rad) + X
    y = 42.0 + t * np.sin(theta_rad) + np.exp(M * np.abs(t)) * np.sin(0.3 * t) * np.cos(theta_rad)
    return x, y

# -------------------------
# Helpers
# -------------------------
def load_data(filepath=DATA_FILENAME):
    if os.path.exists(filepath):
        df = pd.read_csv(filepath)
    else:
        alt = os.path.join("/mnt/data", filepath)
        if os.path.exists(alt):
            df = pd.read_csv(alt)
        else:
            raise FileNotFoundError(f"Data file not found: '{filepath}'")
    if not {'x', 'y'}.issubset(df.columns):
        raise ValueError("CSV must contain 'x' and 'y' columns")
    return df

def heuristic_initial_guess(x, y):
    """Simple linear approximation ignoring oscillatory term to get initial theta and X."""
    Y = (y - 42.0).reshape(-1, 1)
    A = np.hstack([np.ones_like(Y), Y])
    coeffs, *_ = np.linalg.lstsq(A, x, rcond=None)
    a, b = float(coeffs[0]), float(coeffs[1])
    # b ≈ cot(theta) -> theta = atan(1/b)
    theta_rad = atan(1.0 / b)
    theta_deg = degrees(theta_rad)
    X_est = a
    M_est = 0.0
    return theta_deg, theta_rad, M_est, X_est

def sum_min_distance(params, x_points, y_points, inner_maxiter=INNER_MAXITER, xatol=INNER_XATOL):
    """Compute sum of minimal distances from each observed point to the model curve
       by optimizing over t in [T_MIN, T_MAX] for each point.
       params: [theta_deg, M, X]
    """
    theta_deg, M, X = params
    theta_rad = radians(float(theta_deg))
    total = 0.0
    for xi, yi in zip(x_points, y_points):
        def dist_t(t):
            xt, yt = model_xy(t, theta_rad, M, X)
            return float(((xi - xt)**2 + (yi - yt)**2)**0.5)
        res = minimize_scalar(dist_t, bounds=(T_MIN, T_MAX), method='bounded',
                              options={'xatol': xatol, 'maxiter': inner_maxiter})
        total += float(res.fun)
    return float(total)

def project_all(x_obs, y_obs, theta_rad, M, X, inner_maxiter=FINAL_INNER_MAXITER, xatol=FINAL_XATOL):
    n = len(x_obs)
    t_proj = np.empty(n, dtype=float)
    residuals = np.empty(n, dtype=float)
    for i, (xi, yi) in enumerate(zip(x_obs, y_obs)):
        def dist_t(t):
            xt, yt = model_xy(t, theta_rad, M, X)
            return float(((xi - xt)**2 + (yi - yt)**2)**0.5)
        res = minimize_scalar(dist_t, bounds=(T_MIN, T_MAX), method='bounded',
                              options={'xatol': xatol, 'maxiter': inner_maxiter})
        t_proj[i] = float(res.x)
        residuals[i] = float(res.fun)
    return t_proj, residuals

# -------------------------
# Main fitting routine
# -------------------------
def fit_parameters(df):
    x_vals = df['x'].values
    y_vals = df['y'].values
    n = len(df)
    print(f"Data points: {n}")

    # heuristic initialization
    theta_deg0, theta_rad0, M0, X0 = heuristic_initial_guess(x_vals, y_vals)
    print("Initial guess (linear approx):")
    print(f"  theta (deg): {theta_deg0:.6f}, X: {X0:.6f}, M: {M0:.6f}")

    # subsample for DE
    idx = np.linspace(0, n - 1, min(SUBSAMPLE_SIZE, n), dtype=int)
    x_sub = x_vals[idx]
    y_sub = y_vals[idx]

    # prepare objective for DE (use partial to pass subsample)
    de_obj = partial(sum_min_distance, x_points=x_sub, y_points=y_sub,
                     inner_maxiter=INNER_MAXITER, xatol=INNER_XATOL)

    print(f"Running differential evolution on subsample (size={len(x_sub)}) ...")
    start = time.time()
    try:
        de_result = differential_evolution(de_obj, BOUNDS, maxiter=DE_MAXITER, popsize=DE_POPSZ,
                                          tol=1e-6, polish=False, seed=RANDOM_SEED, workers=DE_WORKERS)
    except TypeError:
        # older SciPy may not accept 'workers'
        de_result = differential_evolution(de_obj, BOUNDS, maxiter=DE_MAXITER, popsize=DE_POPSZ,
                                          tol=1e-6, polish=False, seed=RANDOM_SEED)
    dt_de = time.time() - start
    print(f"DE done in {dt_de:.1f} s; best sub-objective = {de_result.fun:.6f}")
    init_params = de_result.x

    # local refinement on full dataset
    print("Refining on full dataset (local optimization)...")
    start = time.time()
    full_obj = partial(sum_min_distance, x_points=x_vals, y_points=y_vals,
                       inner_maxiter=INNER_MAXITER, xatol=INNER_XATOL)
    res_local = minimize(full_obj, init_params, method='L-BFGS-B', bounds=BOUNDS,
                         options={'ftol': 1e-9, 'maxiter': 200})
    dt_local = time.time() - start
    print(f"Local refine done in {dt_local:.1f} s; full objective = {res_local.fun:.6f}")

    theta_deg_fit, M_fit, X_fit = float(res_local.x[0]), float(res_local.x[1]), float(res_local.x[2])
    theta_rad_fit = radians(theta_deg_fit)

    # final projection to compute t_i and residuals with higher inner tolerance
    print("Projecting all points to curve (final pass)...")
    start = time.time()
    t_proj, residuals = project_all(x_vals, y_vals, theta_rad_fit, M_fit, X_fit,
                                    inner_maxiter=FINAL_INNER_MAXITER, xatol=FINAL_XATOL)
    dt_proj = time.time() - start

    sum_L1 = float(residuals.sum())
    mean_r = float(residuals.mean())
    max_r = float(residuals.max())

    # save results
    out_df = pd.DataFrame({
        'x': x_vals,
        'y': y_vals,
        't_proj': t_proj,
        'residual': residuals
    })
    out_df.to_csv(OUT_CSV, index=False)

    results = {
        'theta_deg': theta_deg_fit,
        'theta_rad': theta_rad_fit,
        'M': M_fit,
        'X': X_fit,
        'sum_L1': sum_L1,
        'mean_residual': mean_r,
        'max_residual': max_r,
        'de_time_s': dt_de,
        'local_time_s': dt_local,
        'projection_time_s': dt_proj
    }
    return results

# -------------------------
# Entrypoint
# -------------------------
def main():
    t0 = time.time()
    df = load_data(DATA_FILENAME)
    results = fit_parameters(df)
    total = time.time() - t0

    print("\nFitted parameters:")
    print(f"  theta (deg): {results['theta_deg']:.8f}")
    print(f"  theta (rad): {results['theta_rad']:.8f}")
    print(f"  M:           {results['M']:.8f}")
    print(f"  X:           {results['X']:.8f}")
    print("\nResidual summary:")
    print(f"  sum_L1: {results['sum_L1']:.6f}")
    print(f"  mean:   {results['mean_residual']:.6f}")
    print(f"  max:    {results['max_residual']:.6f}")
    print(f"\nTiming (s): DE={results['de_time_s']:.1f}, local={results['local_time_s']:.1f}, projection={results['projection_time_s']:.1f}")
    print(f"Total elapsed time: {total:.1f} s")
    # produce a short README describing approach
    readme_text = (
        "Submission README\n\n"
        "Method summary:\n"
        "- Linear heuristic used to initialize parameters.\n"
        "- Differential evolution (global search) on an evenly spaced subsample to locate a good region.\n"
        "- Local L-BFGS-B refinement on the full dataset with per-point 1D projection (t in [6,60]).\n\n"
        "Files generated:\n"
        f"- {OUT_CSV}: per-point projected parameter t and residuals.\n\n"
        "How to run:\n"
        "1) Place 'xy_data.csv' in the same folder as this script.\n"
        "2) Run: python fit_curve_submission.py\n\n"
        "Contact: Student submission code (no external generation notes)."
    )
    with open(OUT_README, "w") as fh:
        fh.write(readme_text)
    print(f"\nREADME saved to '{OUT_README}' and results saved to '{OUT_CSV}'.")

if __name__ == "__main__":
    main()
