# Small Wind Turbine Blade-Element Momentum Analysis

## Overview

This project implements a Blade-Element Momentum (BEM) theory solver in MATLAB for the aerodynamic analysis and preliminary design of a small wind turbine rotor.

The model combines blade-element aerodynamics with momentum theory to calculate the performance of a rotor across a range of operating conditions. The implementation includes Prandtl tip and root loss corrections, wake rotation, iterative induction-factor calculations, and the Buhl modification of the Glauert high-induction correction*.

The project also includes tools for generating blade geometry, loading aerofoil polar data, investigating the effect of tip-speed ratio, and comparing different blade designs.

---

## Key Features

- Blade-Element Momentum (BEM) aerodynamic model
- Iterative solution for axial and tangential induction factors
- Prandtl tip-loss correction
- Prandtl root-loss correction
- Wake rotation through tangential induction
- Buhl modification of the Glauert high-induction correction
- Aerofoil polar data interpolation
- Parametric blade geometry generation
- Tip-speed-ratio sweep
- Power and thrust coefficient calculations
- Starting-torque analysis
- Blade-design sensitivity analysis

---

## Methodology

The rotor is divided into a number of radial blade elements. For each element, the local flow conditions are determined from the axial and tangential induction factors.

The solver iteratively calculates:

1. Local inflow angle
2. Local angle of attack
3. Aerofoil lift and drag coefficients
4. Normal and tangential force coefficients
5. Prandtl tip and root loss factors
6. Axial induction factor, `a`
7. Tangential induction factor, `a'`

The induction factors are updated iteratively using a relaxation factor until the solution converges.

The resulting elemental thrust and torque contributions are then integrated across the blade to obtain the overall rotor performance.

The core solver returns quantities including:

- Axial induction factor, `a`
- Tangential induction factor, `a'`
- Inflow angle, `φ`
- Local angle of attack, `α`
- Lift coefficient, `CL`
- Drag coefficient, `CD`
- Prandtl loss factor, `F`
- Elemental thrust, `dT`
- Elemental torque, `dQ`
- Total thrust, `T`
- Total torque, `Q`
- Power, `P`
- Power coefficient, `CP`
- Thrust coefficient, `CT`
- Tip-speed ratio, `λ`

---

## Blade Geometry

The blade geometry is generated parametrically from a set of design parameters:

- Rotor radius, `R`
- Root radius, `r0`
- Number of blades, `B`
- Design tip-speed ratio, `λ`
- Design angle of attack, `α_des`
- Root-to-tip chord ratio
- Number of blade elements

The blade planform is generated using a linear chord distribution, while the aerodynamic twist distribution is calculated from the specified design tip-speed ratio and design angle of attack.

The geometry generation routine calculates the blade planform area, root chord, tip chord, and twist distribution.

---

## Aerofoil Polar Data

The solver uses aerofoil polar data to obtain the sectional lift and drag coefficients as a function of angle of attack.

The `loadPolar.m` function reads CSV files containing:

- Angle of attack, `α`
- Lift coefficient, `CL`
- Drag coefficient, `CD`

The data is sorted, cleaned, and converted into interpolation functions for use by the BEM solver.

The current example uses an SG6043 aerofoil polar at approximately Re = 100,000. The polar filename is specified in `main_BEM.m`.

---

## Rotor Configuration

The current example uses the following design parameters:

| Parameter | Value |
|---|---:|
| Rotor radius, `R` | 0.25 m |
| Root radius, `r0` | 0.05 m |
| Number of blades, `B` | 2 |
| Design tip-speed ratio, `λ` | 3 |
| Design angle of attack | 6° |
| Root/tip chord ratio | 2.5 |
| Blade elements | 20 |
| Design wind speed | 10 m/s |
| Air density | 1.225 kg/m³ |

These parameters can be modified directly in `main_BEM.m`.

---

## Analysis

### Design-Point Analysis

The main script first generates the blade geometry and evaluates the rotor at the selected design operating point.

The analysis produces information including:

- Blade chord distribution
- Blade twist distribution
- Axial induction distribution
- Tangential induction distribution
- Local angle of attack
- Aerofoil lift and drag coefficients
- Prandtl loss factor
- Power coefficient, `CP`
- Thrust coefficient, `CT`
- Rotor power
- Rotor torque

The implementation also reports the number of blade elements that successfully converged. 

### Tip-Speed-Ratio Sweep

The `sweepLambda.m` function evaluates the rotor over a range of tip-speed ratios while maintaining a constant wind speed.

This allows the rotor's power coefficient, `CP`, to be plotted against tip-speed ratio and the operating point corresponding to maximum power coefficient to be identified. 

The current implementation evaluates:

λ = 0.5 → 8
