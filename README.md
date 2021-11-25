# lbm

![master badge](https://github.com/jviquerat/lbm/workflows/lbm/badge.svg?branch=master)

A simple lattice-Boltzmann code for 2D flow resolutions. All the tools are contained in the `lattice.py` file, and separate cases are built on top of this library.

<p align="center">
  <img width="900" alt="" src="https://user-images.githubusercontent.com/44053700/99295075-3dd29f00-2845-11eb-8e05-1d8a132b0feb.gif">
</p>

## Contents

This LBM code includes:

- D2Q9 lattice
- TRT collision operator
- Zou-He on all boundary conditions
- Drag/lift computation using interpolated bounce-back
- Core routines are deferred to Numba

## Running simulations

Cases are described in the `lbm/src/app/` repository. To run a simulation, adjust the parameters in the related python file, then run `python3 start.py <app_name>`. A results folder will be generated in `./results/` with the current date and time. If you wish to add a new application, you must create a new app, and register it in the factory, located in `lbm/src/app/app.py`. Below are some examples and benchmarks that were ran with the code. The related cases are available in the repository.

## Benchmarks

### Lid-driven cavity

A simple driven cavity in unit square. Below are the computed time-domain velocity norms and final streamlines at Re=100 (left) and Re=1000 (right).

<p align="center">
  <img width="350" alt="" src="lbm/save/driven_cavity/re_100_nx_200/anim-opt.gif"> <img width="350" alt="" src="lbm/save/driven_cavity/re_1000_nx_250/anim-opt.gif">
</p>

A comparison of `u = f(y)` and `v = f(x)` at the center of the domain with reference data from <a href="https://www.sciencedirect.com/science/article/pii/0021999182900584">"U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow using Navier-Stokes equations and multigrid method*"</a>.

<p align="center">
  <img width="350" alt="" src="lbm/save/driven_cavity/re_100_nx_200/re_100.png"> <img width="350" alt="" src="lbm/save/driven_cavity/re_1000_nx_250/re_1000.png">
</p>

### Turek benchmark

The Turek cylinder benchmark CFD case is described in <a href="https://link.springer.com/chapter/10.1007/978-3-322-89849-4_39">"Schafer, M., Turek, S. *Benchmark Computations of Laminar Flow Around a Cylinder*"</a>. The 2D case consists in a circular cylinder in a channel with top and bottom no-slip conditions, and with a Poiseuille flow at the inlet (these cases are named 2D-1 and 2D-2 in the aforementionned reference). The cylinder is voluntarily not centered to trigger instability at sufficient Reynolds number. Here, we explore the accuracy of the drag and lift computation (using IBB). Note that for the 2D-2 case, the values correspond to the max drag and lift. The computational times are obtained on a standard laptop, with a field output every 500 iterations.

<center>

|        |`ny` | 2D-1 (Re=20) Cd, Cl, CPU   | 2D-2 (Re=100) Cd, Cl, CPU  |
|--------|-----|----------------------------|----------------------------|
| Turek  | --- |  5.5800 - 0.0107 - N/A     |  3.2300 - 1.0000 - N/A     |
| lbm    | 100 |  6.6285 - 0.1179 - 238 s.  |  2.8573 - 1.8082 - 384 s.  |
| lbm    | 200 |  5.6139 - 0.0173 - 1431 s. |  3.2247 - 1.1176 - 2493 s. |

</center>

Below is a video of the 2D-2 case:

<p align="center">
  <img width="800" alt="" src="lbm/save/turek/re_100_ny_200/turek.gif">
</p>

## Applications

### Poiseuille with random obstacles

It is possible to run a Poiseuille flow with random obstacles in the domain. Below is an example.

<p align="center">
  <img width="900" alt="" src="https://user-images.githubusercontent.com/44053700/99222254-768d5c80-27e2-11eb-96a0-c26ecadfa0c0.gif">
</p>

### Array of obstacles
