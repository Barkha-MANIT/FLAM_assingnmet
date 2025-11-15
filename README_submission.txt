Submission README

Method summary:
- Linear heuristic used to initialize parameters.
- Differential evolution (global search) on an evenly spaced subsample to locate a good region.
- Local L-BFGS-B refinement on the full dataset with per-point 1D projection (t in [6,60]).

Files generated:
- fit_results.csv: per-point projected parameter t and residuals.

How to run:
1) Place 'xy_data.csv' in the same folder as this script.
2) Run: python fit_curve_submission.py

Contact: Barkha Verma, devbarkha07@gmail.com.