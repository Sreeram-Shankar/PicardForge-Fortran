# 📦 PicardForge-Fortran
### *A fixed-step ODE & PDE time-integration library featuring Picard–Gauss–Seidel implicit solvers and classical explicit methods — implemented in modern Fortran (F2008+)*

---

## ✨ Overview

**PicardForge-Fortran** is a high-performance fixed-step ODE solver library written entirely in modern Fortran.  
It mirrors the design of PicardForge (Python) and PicardForge-Julia, but is optimized for speed, static typing, and efficient array-level operations.

It contains a full suite of **explicit and implicit integrators**, all designed for **semi-discretized PDEs** such as:

- heat/diffusion equations  
- conduction models  
- parabolic finite-difference / finite-volume systems  
- large stiff ODE systems arising from PDE discretization  

All implicit solvers use **Picard fixed-point iteration with Gauss–Seidel relaxation**, eliminating the need for Jacobians or Newton iterations.

This makes PicardForge-Fortran extremely lightweight and fast for diffusion-dominated PDEs.

---

## 🚀 Features

### ✔ Comprehensive ODE solver suite

| Family | Methods | Notes |
|--------|---------|-------|
| **Explicit RK** | RK1–RK5 | Classic explicit Butcher tables |
| **Adams–Bashforth** | AB2–AB5 | Explicit multistep |
| **Adams–Moulton** | AM2–AM5 | Implicit multistep (Picard–GS) |
| **BDF** | BDF1–BDF6 | Fully implicit, robust for stiff PDEs |
| **SDIRK** | SDIRK2–SDIRK4 | Diagonally implicit RK |
| **IRK Collocation** | Gauss, Radau IIA, Lobatto IIIC (s=1–5) | Fully implicit, A- or L-stable |

All implemented with:

- modern Fortran modules  
- assumed-shape arrays  
- explicit interfaces  
- pure matrix-free Picard iteration  

---

## 📁 Repository Structure

```

PicardForge-Fortran/
│
├── ab_module.f90        # AB2–AB5 multistep
├── am_module.f90        # AM2–AM5 implicit multistep (Picard–GS)
├── bdf_module.f90       # BDF1–BDF6 implicit (Picard–GS)
├── rk_module.f90        # RK1–RK5 explicit methods
├── sdirk_module.f90     # SDIRK2–SDIRK4
├── irk_module.f90       # Gauss/Radau/Lobatto collocation IRK
└── picardforge.f90      # Unified module exporting all solvers

````

All solver modules follow the same API style:

```fortran
subroutine solve_method(f, t0, tf, y0, h, n, Yout, tgrid, [args])
````

Where:

* `f` is a user-supplied ODE RHS subroutine
* `y0(:)` is the initial state
* `Yout(k,:)` contains the solution at step k
* fixed-step time integration is used

---

## 🧠 How Picard–Gauss–Seidel Nonlinear Iteration Works

All implicit methods reduce to solving:

[
Y = y_n + h A F(Y)
]

where:

* ( A ) is the Butcher matrix (IRK, SDIRK)
* or derived coefficients (BDF, Adams–Moulton)
* ( F(Y) ) computes stage derivatives

Instead of Newton iteration, PicardForge-Fortran uses:

### **Picard fixed-point iteration**

[
Y^{(k+1)} = y_n + h A F(Y^{(k)})
]

### **With Gauss–Seidel stage relaxation**

Stages update one at a time:

```fortran
do sweep = 1, sweeps
    Y_old = Y_stages
    do i = 1, s
        rhs = 0.0d0
        do j = 1, s
            call f(t + c(j)*h, Y_stages(j,:), fval)
            rhs = rhs + A(i,j) * fval
        end do
        Y_stages(i,:) = y + h * rhs
    end do
    if (norm2(Y_stages - Y_old) < tol) exit
end do
```

### Why this works so well

* Parabolic PDEs have smoothing operators
* Picard iteration converges rapidly with these operators
* Gauss–Seidel makes the iteration more stable
* No Jacobians, no matrix solves, no GMRES — just pure local evaluation

Perfect for large PDE systems from FD/FV/FEM discretization.

---

## 🔧 Usage Examples

### 1️⃣ Solving a simple ODE (exponential decay)

```fortran
program test_rk4
    use rk_module
    implicit none

    real(8) :: y0(1), t0, tf, h
    real(8), allocatable :: Yout(:,:), tgrid(:)

    y0 = [1.0d0]
    t0 = 0.0d0
    tf = 10.0d0
    h  = 0.1d0

    allocate(Yout(int((tf-t0)/h)+1, 1))
    allocate(tgrid(int((tf-t0)/h)+1))

    call solve_rk4(f, t0, tf, y0, h, 1, Yout, tgrid)

contains
    subroutine f(t, y, dydt)
        real(8), intent(in) :: t, y(:)
        real(8), intent(out) :: dydt(:)
        dydt(1) = -0.01d0 * y(1)
    end subroutine f
end program test_rk4
```

---

### 2️⃣ Semi-discretized PDE (heat equation)

```fortran
call solve_bdf4(rhs_heat, 0.0d0, 1.0d0, T0, h, N, Tout, tgrid, sweeps=20, tol=1e-10)
```

Where `rhs_heat(t, T, dTdt)` applies your discrete Laplacian.

---

### 3️⃣ IRK collocation (Radau IIA)

```fortran
call solve_collocation(f, 0.0d0, 1.0d0, y0, 0.05d0, n, Yout, tgrid, "radau", 3, 20, 1e-10)
```

Radau IIA is L-stable → extremely good for stiff diffusion PDEs.

---

## 📈 Stability Notes

### Explicit methods

* **RK and AB** obey CFL-like limits
* Not suitable for stiff PDEs

### Implicit methods

* **AM2** is A-stable
* **BDF family** is excellent for stiff systems (BDF2–BDF6 widely used in PDE solvers)
* **Gauss IRK** → symplectic, high accuracy
* **Radau IIA** → **L-stable**, best for very stiff PDEs
* **Lobatto IIIC** → symmetric and stiffly accurate

### Picard–GS

* Fast convergence for diffusion problems
* No matrix solves
* Zero Jacobian work

---

## 🛠 Building & Requirements

* Fortran **2008+**
* A modern compiler: `gfortran`, `ifort`, `nvfortran`, etc.
* No external dependencies

Example build:

```bash
gfortran -O3 -c rk_module.f90 am_module.f90 ab_module.f90 bdf_module.f90 irk_module.f90 sdirk_module.f90 picardforge.f90
```

---

## 📜 License

Released under the **MIT License**.

---

## 🤝 Contributing

Contributions welcome! Ideas include:

* Adaptive time-stepping
* Embedded IRK error estimators
* Multi-rate or IMEX schemes
* GPU offloading (OpenMP target / CUDA Fortran)
* Parallel domain decomposition for PDEs

Enjoy high-performance implicit/explicit integration with **PicardForge-Fortran** 🚀

```
Just tell me.
```
