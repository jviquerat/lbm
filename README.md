# lbm
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

Below are some examples ran with the code. The related cases are available in the repository.

## Lid-driven cavity

A simple driven cavity in unit square. Launch it by running ```python3 cavity.py```.  
Below are the computed time-domain velocity norms and final streamlines at Re=100 (left) and Re=1000 (right).

<p align="center">
  <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/98362755-12260c80-202e-11eb-873a-7f4948150331.gif"> <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/98362762-1520fd00-202e-11eb-8495-b21ecac968bc.gif">
</p>

<p align="center">
  <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/76288545-4088f780-62a7-11ea-9893-dd0a19339bc5.png"> <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/76447230-9a351300-63c8-11ea-8722-35e1eb2151c0.png">
</p>

A comparison of u_x = f(y) at the center of the domain with reference data from "U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow using Navier-Stokes equations and multigrid method*."

<p align="center">
  <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/76288543-3ebf3400-62a7-11ea-9e2b-13e0f6327c89.png"> <img width="350" alt="" src="https://user-images.githubusercontent.com/44053700/76447238-9e613080-63c8-11ea-8e60-6f77248518a2.png">
</p>

## Turek benchmark

The Turek cylinder benchmark CFD case is described in "Schafer, M., Turek, S. *Benchmark Computations of Laminar Flow Around a Cylinder*". The 2D case consists in a circular cylinder in a channel with top and bottom no-slip conditions, and with a Poiseuille flow at the inlet (these cases are named 2D-1 and 2D-2 in the aforementionned reference). The cylinder is voluntarily not centered to trigger instability at sufficient Reynolds number. Here, we explore the accuracy of the drag and lift computation.

|        |`ny` | 2D-1 (Re=20) Cd   | 2D-1 (Re=20) Cl   | 2D-2 (Re=100) Cd  | 2D-2 (Re=100) Cl  |
|--------|-----|-------------------|-------------------|-------------------|-------------------|
| Turek  | --- |  5.5800           |  0.0107           |  3.2300           |  1.0000           |
| lbm    | 100 |  5.6300           |  0.0862           |  3.0411           |  0.5834           |
| lbm    | 200 |  5.5804           |  0.0371           |  3.2582           |  1.2047           |
| lbm    | 300 |  5.5846           |  0.0261           |  3.2152           |  1.0987           |

Below are videos of the 2D-1 and 2D-2 cases:

<p align="center">
  <img width="800" alt="" src="https://user-images.githubusercontent.com/44053700/101684500-d0421900-3a66-11eb-8137-7d936569c388.gif">
  <img width="800" alt="" src="https://user-images.githubusercontent.com/44053700/101684599-f7004f80-3a66-11eb-820d-c4da29dc1dfe.gif">
</p>

## Poiseuille with random obstacles

It is possible to run a Poiseuille flow with random obstacles in the domain. Below is an example.

<p align="center">
  <img width="900" alt="" src="https://user-images.githubusercontent.com/44053700/99222254-768d5c80-27e2-11eb-96a0-c26ecadfa0c0.gif">
</p>

## Running

To run a simulation, adjust the parameters in the related python file, then run ```python3 case.py```. A results folder will be generated in ```./results/``` with the current date and time. The ```png/``` folder will contain outputs of the velocity norm over the domain. To generate a video out of the png files, you can use the ```convert``` command as follows:

```convert -delay 10 -resize 50% -loop 0 'u_norm_%d.png'[0-100] animation.gif```

To optimize and resize gifs, use ```gifsicle``` :

```gifsicle -i animation.gif --scale 0.6 -O3 --colors 256 -o anim-opt.gif```
