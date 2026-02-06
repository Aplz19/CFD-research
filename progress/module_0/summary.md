# Module 0 Summary

**Status:** Completed

**Purpose:** Build reusable infrastructure for all later CFD modules

---

## Key Components

### 1. SSPRK Time Integrators
Solve semi-discrete systems of the form:
```
dU/dt = L(U, t)
```
where L is the spatially discretized residual/RHS operator.

**Forward Euler step operator:**
```
F(U, t, dt) = U + dt * L(U, t)
```

**SSPRK2 (2nd order, 2 stages):**
```
U^(1) = F(U^n, t^n, dt)
U^{n+1} = 0.5*U^n + 0.5*F(U^(1), t^n + dt, dt)
```

**SSPRK3 (3rd order, 3 stages):**
```
U^(1) = F(U^n, t^n, dt)
U^(2) = 0.75*U^n + 0.25*F(U^(1), t^n + dt, dt)
U^{n+1} = (1/3)*U^n + (2/3)*F(U^(2), t^n + 0.5*dt, dt)
```

### 2. Grid Files (Structured Nodes)

**1D Format:**
- Line 1: Nx (number of nodes)
- Lines 2 to Nx+1: x(i) coordinates
- Cell i spans [x_i, x_{i+1}], so Nx nodes = Nx-1 cells

**2D Format:**
- Line 1: Nx Ny
- Next Nx*Ny lines: x(i,j) y(i,j) with j outer loop, i inner loop
- (Nx+1) x (Ny+1) nodes = Nx x Ny cells

### 3. Geometry Preprocessing

**1D:** From node coordinates x_{i+1/2}, compute:
- Cell centers: x_i = 0.5*(x_{i+1/2} + x_{i-1/2})
- Cell widths: dx_i = x_{i+1/2} - x_{i-1/2}

**2D:** From corner nodes, compute:
- Cell centers (xc, yc)
- Cell areas using shoelace formula
- Outward unit normals at each face

---

## Required MATLAB Interfaces

**Time Steppers:**
- `U = ssprk2(U, t, dt, rhs, geom, params)`
- `U = ssprk3(U, t, dt, rhs, geom, params)`

**Grid I/O:**
- `grid = read_grid_1d(filename)`
- `grid = read_grid_2d(filename)`

**Geometry:**
- `geom = build_geom_1d(grid)` returns xc, dx
- `geom = build_geom_2d(grid)` returns xc, yc, A, face vectors

**Residual signature:**
```
R = rhs(U, t, geom, params)
```

---

## Examples Completed

**Example A: Nonautonomous ODE**
- Equation: u'(t) = -u + sin(t), u(0) = 0, t in [0, 10]
- Exact solution: u(t) = 0.5*(sin(t) - cos(t)) + 0.5*e^{-t}
- Tasks: Test SSPRK2/SSPRK3 with multiple dt, compare to ode45, report errors

**Example B: 2D Grid Preprocessing**
- Read 2D structured grid
- Compute cell centers, areas, face normals
- Visualize grid with cell centers and face normals overlaid
