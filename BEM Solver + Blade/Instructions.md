# Wind turbine BEM solver (MATLAB)

Full blade-element momentum theory solver for the small wind turbine
project. Includes Prandtl tip + root loss, wake rotation (a'), and the
Buhl modification of Glauert's heavily-loaded rotor correction.

## File layout

| File | Purpose |
| --- | --- |
| `main_BEM.m` | Driver script. Edit the parameters block at the top, then run. |
| `buildBladeGeom.m` | Converts (R, lambda, B, alpha_des, taper) into chord and twist tables, using equations 29 and 31 from the project brief. |
| `loadPolar.m` | Reads aerofoil polar CSVs (AirfoilTools format auto-detected). Returns interpolants for C_L(alpha) and C_D(alpha). |
| `bemSolve.m` | Core solver. Iterates a and a' for each blade element. Returns per-element diagnostics plus integrated C_P, C_T, P, Q. |
| `sweepLambda.m` | Repeats bemSolve over a range of tip speed ratios at fixed wind speed. Used to find the actual peak-C_P operating point. |

## How to use

1. Get an aerofoil polar. Go to airfoiltools.com, find your aerofoil
   (e.g. SG6043), pick a Reynolds number close to your operating range
   (start with Re = 100,000), and download the CSV. Place it in this
   folder.

2. Edit the parameters block at the top of `main_BEM.m`. The starting
   values match the design in the project notes: R = 0.25 m, B = 3,
   lambda = 3, alpha_des = 6 deg, taper = 2.5.

3. Run `main_BEM`. You should get:
   - Geometry plots (chord and twist vs r)
   - Per-element diagnostics at the design point
   - A C_P–lambda curve with the peak marked
   - A table of starting torques at low Omega

## What to look for

- **C_P at design point** should be in the range 0.35–0.45 for a small
  3-blade rotor. Above 0.5 is suspicious (likely a bug or missing tip
  loss). Below 0.25 means your aerofoil is mis-matched or your alpha_des
  is wrong.
- **C_P-lambda peak** may not occur at the lambda you designed for —
  that's fine and useful. If the peak is far from your design lambda,
  consider redesigning the blade for the lambda where the peak is.
- **Starting torque** at Omega -> 0 must exceed 0.03 Nm. If not, either
  increase blade area (lower lambda design point) or use an aerofoil
  with better low-Re post-stall behaviour.
- **Convergence**: per-element `converged` should be true everywhere.
  If not, increase `maxIt` or lower `relax`.

## Caveats and assumptions

- The polar must cover the angle-of-attack range you actually hit.
  At very low Omega (starting) alpha can exceed 30 deg — XFOIL polars
  rarely go that high, and BEM near stall is unreliable anyway. Treat
  starting-torque numbers as indicative, not predictive.
- Only one polar Re is used per run. In reality Re varies along the
  blade. For a more accurate run, generate polars at several Re values
  and interpolate by radius — this is a worthwhile upgrade for the
  final report.
- Glauert/Buhl correction kicks in only when a > 0.4. If you see many
  elements with a near 0.95 (the clamp value), something is wrong —
  almost always a stalled element with negative C_L lookup.
