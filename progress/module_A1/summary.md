# Module A1 Summary

**Status:** Completed

**Topic:** 1st-order Finite Volume method for 1D linear advection equation

---

## Governing Equation

**Linear advection:**
```
u_t + a * u_x = 0
```
where `a` is constant advection speed.

**Exact solution:** Translation of initial condition
```
u(x, t) = u0(x - a*t)
```
With periodic BCs on [0, L], interpret (x - a*t) mod L.

---

## Finite Volume Formulation

### Grid Indexing
- Cell i spans [x_{i-1/2}, x_{i+1/2}]
- Cell center: x_i
- Cell width: dx_i = x_{i+1/2} - x_{i-1/2}
- Supports both uniform and nonuniform grids

### Semi-Discrete Form
```
d(u_bar_i)/dt = -(1/dx_i) * (f_hat_{i+1/2} - f_hat_{i-1/2})
```
where f_hat are numerical fluxes at faces.

### Numerical Flux
Physical flux: f(u) = a*u

**First-order face states:**
- u_L = u_bar_i
- u_R = u_bar_{i+1}

**Local Lax-Friedrichs (LLF) flux:**
```
f_hat(u_L, u_R) = 0.5*(f(u_L) + f(u_R)) - 0.5*alpha*(u_R - u_L)
```
For linear advection with constant a: alpha = |a|

---

## Time Integration

Use SSPRK2 or SSPRK3 from Module 0.

**CFL condition:**
```
dt = CFL * min_i(dx_i) / |a|
```
Typical CFL values: 0.2 to 0.9

---

## Boundary Conditions

**Periodic:**
```
u_bar_0 = u_bar_Nx    (left ghost = rightmost cell)
u_bar_{Nx+1} = u_bar_1  (right ghost = leftmost cell)
```

Boundary fluxes:
```
f_hat_{1/2} = f_hat(u_bar_Nx, u_bar_1)
f_hat_{Nx+1/2} = f_hat(u_bar_Nx, u_bar_1)
```

---

## Required MATLAB Interface

```matlab
R = rhs_advection_1d(U, t, geom, params)
```
- geom contains dx (cell widths), Nx, Nf, xc, xf
- params contains a (advection speed) and CFL

## Implementation Progress

### Files
- `READGRID_1D.m` — reads 1D grid file (Nf, xf)
- `GEOM_1D.m` — computes cell widths (dx), cell centers (xc), stores Nf, Nx, xf
- `rhs_advection_1d.m` — RHS function: physical flux, LLF numerical flux with periodic BCs, residual computation
- `ssprk2.m` — SSPRK2 time integrator using function handles
- `EXA1.m` — main driver script: grid setup, IC, time loop, error computation, plotting for all cases

### Results
- Error table generated for smooth uniform case across CFL = 0.5, 0.9, 1.5 and T = 0.5, 1.0, 1.5
- First-order convergence confirmed: error halves as Nx doubles
- CFL = 0.5 and 0.9 give nearly identical errors; CFL = 1.5 is unstable (blows up)
- Non-smooth case shows numerical diffusion smearing the discontinuity; narrows with more cells but never truly sharp
- Report completed: `EXA1_Cardozo.docx`/`.pdf`

---

## Test Cases

**Domain:** x in [0, 1], a = 1, periodic BCs

**Case 1 - Smooth IC:**
```
u0(x) = sin(2*pi*x)
```
- Grid: Nx = 2^{i+3} for i = 1, 2, ..., 9 (i.e., 16, 32, 64, ..., 4096 cells)
- Both uniform and nonuniform grids

**Case 2 - Discontinuous (top-hat pulse):**
```
u0(x) = 1  if 0.25 <= x <= 0.5
        0  otherwise
```
- Grid: Nx = 100, 500, 1000, 2500, 5000
- Both uniform and nonuniform grids

**Final times:** T in [0.25, 1.5], e.g., T = 0.5 (half period), T = 1 (full period)

---

## Report Requirements

1. **Plots:** Overlay u(x,0), u(x,T), u_exact(x,T) at cell centers

2. **Error table (smooth case, uniform grid):**
   ```
   ||e||_2 = sqrt(sum_i (u_bar_i(T) - u_exact(x_i, T))^2 * dx_i)
   ```

3. **Discussion points:**
   - How FV captures smooth solutions as Nx increases
   - How FV captures discontinuities as Nx increases
   - What happens when CFL = 1.5 (instability demonstration)

---

## Key MATLAB Snippets

**Grid setup (uniform/nonuniform):**
```matlab
if uniform
    xf = linspace(0, L, Nx+1).';  % face coordinates
else
    s = linspace(0, 1, Nx+1).';
    beta = 2.0;  % stretching
    xf = L*(exp(beta*s) - 1)/(exp(beta) - 1);
end
dx = diff(xf);
xc = 0.5*(xf(1:end-1) + xf(2:end));  % cell centers
```

**CFL timestep:**
```matlab
dt = CFL * min(dx) / abs(a);
```
