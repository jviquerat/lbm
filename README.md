# lbm
A simple lattice-Boltzmann code for 2D channel flows

![animation](https://user-images.githubusercontent.com/44053700/73072734-0bbb0f80-3eb6-11ea-8703-15145838ede7.gif)

## Contents

This LBM code uses some of the most basic assumptions:

- D2Q9 lattice
- BGK collision operator
- Zou-He on all boundary conditions, i.e. prescribed $$u_x$$ and uy at inlet, prescribed rho and uy at outlet, prescribed ux and uy on top and bottom walls (no-slip condition)

As of now, it is limited to channel flows with an obstacle. The ```Shape``` is used to generate random shapes (or cylinder, or any shape that can be read from an in-house ```.csv``` format, see here https://github.com/jviquerat/bezier_shapes). The lattice is then generated, the given shape being centered on ```(0,0)```.

## Running

To run a simulation, adjust the parameters in ```start.py```, then run ```python3 start.py```. A results folder will be generated in ```./results/``` with the current date and time. The ```png/``` folder will contain outputs of the velocity norm over the domain. To generate a video out of the png files, you can use the ```convert``` command as follows:

```convert -delay 10 -loop 0 'vel_%d.png'[0-100] animation.gif```
