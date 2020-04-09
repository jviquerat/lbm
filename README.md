# lbm
A simple lattice-Boltzmann code for 2D flows resolutions

<img width="900" alt="" src="https://user-images.githubusercontent.com/44053700/78360815-a8391680-75b7-11ea-9a4b-98114fbf8bb9.gif">

## Contents

This LBM code includes:

- D2Q9 lattice
- TRT collision operator
- Zou-He on all boundary conditions
- Drag/lift computation using interpolated bounce-back

Below are some examples ran with the code. The related cases are available in the repository.

## Lid-driven cavity

A simple driven cavity in unit square. Launch it by running ```python3 cavity.py```

Below are the computed time-domain velocity norms at Re=100 (left) and Re=1000 (right).

<p align="center">
<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/77251271-bac55e80-6c4d-11ea-9c75-bdc10da0fef9.gif"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447194-8ab5ca00-63c8-11ea-8d80-0fc7f9c17ed2.gif">
</p>

The final streamlines:

<p align="center">
<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76288545-4088f780-62a7-11ea-9893-dd0a19339bc5.png"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447230-9a351300-63c8-11ea-8722-35e1eb2151c0.png">
</p>

A comparison of u_x = f(y) at the center of the domain with reference data from Ghia

<p align="center">
<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76288543-3ebf3400-62a7-11ea-9e2b-13e0f6327c89.png"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447238-9e613080-63c8-11ea-8e60-6f77248518a2.png">
</p>

## Poiseuille flow

The establishment of a Poiseuille flow in a channel. Launch it by running ```python3 poiseuille.py```

First, the evolution of velocity norm with time :

<p align="center">
<img width="550" alt="" src="https://user-images.githubusercontent.com/44053700/77248108-4b447480-6c37-11ea-8396-31207aad9bc8.gif">
</p>

And the comparison of u_x = f(y) at the center of the domain with exact solution:

<p align="center">
<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/77248104-47185700-6c37-11ea-8e2d-693e34a0132c.png">
</p>

## Poiseuille with obstacle

The computation of drag and lift force on an arbitrary shape in a Poiseuille flow. Launch it by running ```python3 poiseuille.py```

<img width="900" alt="" src="https://user-images.githubusercontent.com/44053700/78360815-a8391680-75b7-11ea-9a4b-98114fbf8bb9.gif">

## Running

To run a simulation, adjust the parameters in the related python file, then run ```python3 case.py```. A results folder will be generated in ```./results/``` with the current date and time. The ```png/``` folder will contain outputs of the velocity norm over the domain. To generate a video out of the png files, you can use the ```convert``` command as follows:

```convert -delay 10 -loop 0 'u_norm_%d.png'[0-100] animation.gif```

To optimize and resize gifs, use ```gifsicle``` :

```gifsicle -i animation.gif --scale 0.6 -O3 --colors 256 -o anim-opt.gif```
