# lbm
A simple lattice-Boltzmann code for simple 2D channel applications

<!---- ![animation](https://user-images.githubusercontent.com/44053700/73072734-0bbb0f80-3eb6-11ea-8703-15145838ede7.gif) -->

## Contents

This LBM code includes:

- D2Q9 lattice
- TRT collision operator
- Zou-He on all boundary conditions

Below are some examples ran with the code. The related cases are available in the repository.

## Lid-driven cavity

A simple driven cavity in unit square. Launch it by running ```python3 cavity.py```

*Top*: Re=100  
*Bottom*: Re=1000  
*Left*: evolution of velocity norm with time.   
*Center*: velocity streamlines at convergence.   
*Right*: comparison of u_x = f(y) at the center of the domain with reference data from Ghia

<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/77251271-bac55e80-6c4d-11ea-9c75-bdc10da0fef9.gif"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76288545-4088f780-62a7-11ea-9893-dd0a19339bc5.png"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76288543-3ebf3400-62a7-11ea-9e2b-13e0f6327c89.png">

<img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447194-8ab5ca00-63c8-11ea-8d80-0fc7f9c17ed2.gif"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447230-9a351300-63c8-11ea-8722-35e1eb2151c0.png"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/76447238-9e613080-63c8-11ea-8e60-6f77248518a2.png">

## Poiseuille flow

The establishment of a Poiseuille flow in a channel at Re=100. Launch it by running ```python3 poiseuille.py```

<img width="550" alt="" src="https://user-images.githubusercontent.com/44053700/77248108-4b447480-6c37-11ea-8396-31207aad9bc8.gif"> <img width="280" alt="" src="https://user-images.githubusercontent.com/44053700/77248104-47185700-6c37-11ea-8e2d-693e34a0132c.png">

## Poiseuille with obstacle

<!--- As of now, it is limited to channel flows with an obstacle. The ```Shape``` is used to generate random shapes (or cylinder, or any shape that can be read from an in-house ```.csv``` format, see here https://github.com/jviquerat/bezier_shapes). The lattice is then generated, the given shape being centered on ```(0,0)```. -->

## Running

To run a simulation, adjust the parameters in ```start.py```, then run ```python3 start.py```. A results folder will be generated in ```./results/``` with the current date and time. The ```png/``` folder will contain outputs of the velocity norm over the domain. To generate a video out of the png files, you can use the ```convert``` command as follows:

```convert -delay 10 -loop 0 'u_norm_%d.png'[0-100] animation.gif```
