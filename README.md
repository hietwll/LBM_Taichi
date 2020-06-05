# LBM_Taichi

This script implements a 2d fluid solver based on [Lattice Boltzmann method](https://en.wikipedia.org/wiki/Lattice_Boltzmann_methods) using [Taichi](https://github.com/taichi-dev/taichi) programming language. The high-performance cross-platform CFD (computational fluid dynamics) solver can be achieved within 200 lines thanks to taichi.


## Usage
To numerically solve a fluid-dynamics problem, the domain size, fluid property, boundary conditions and initial conditions should be given. In this code, these parameters can be specified by instancing the solver:
```
lbm = lbm_solver(nx, ny, niu, bc_type, bc_value)
```
The meanning of each parameter is:
- ``nx``, ``ny`` are domain size. Note they are given in dimensionless form (ie. lattice units), which assumes ``dx = dy = dt = 1.0``, where ``dx`` and ``dy`` are discretized grid sizes, ``dt`` is the time interval of one step.
- ``niu`` is the fluid viscosity in lattice units. Note there is a transformation between SI units and lattice units.
- ``bc_type`` is a four-element python list denoting the ``[left, top, right, bottom]`` boundary condition type. The velocity at the boundary is set based on ``bc_type``. If ``bc_type = 0``, velocity is set as constant value (Dirichlet condition ) given in ``bc_value``. If ``bc_type = 1``, the derivative of velocity in boudary normal direction is set to zero (Neumann condition).
- ``bc_value`` is a ``(4,2)`` python list, it gives constant velocity value at each boudary when ``bc_type = 0``.

## Example1: Lid-driven Cavity Flow
<div align="center">
<img src="https://raw.githubusercontent.com/hietwll/common_files/master/graphics/lbm_taichi/LidDrivenCavity.png" height="400px">
</div>

Lid-driven cavity flow is benchmark fluid-dynamics problem used to verificate the solver. To compare simulation results based on different unit-systems, the flow Reynolds number ``Re`` should keep the same. For this flow, ``Re`` is defined as ``Re = U * L / niu``, so a solver with `` Re = 1000 `` can be given by:
```
lbm = lbm_solver(256, 256, 0.0255, [0, 0, 0, 0], 
      [[0.0, 0.0], [0.1, 0.0], [0.0, 0.0], [0.0, 0.0]])
```
Here ``Re = U * (nx-1) * dx / niu = 0.1 * 255.0 / 0.0255``. The velocity magnitude is shown in the contour below and x-componet of velocity in the middle line is compared with literature results.

<img src="https://raw.githubusercontent.com/hietwll/common_files/master/graphics/lbm_taichi/lid.gif" height="293px"> <img src="https://raw.githubusercontent.com/hietwll/common_files/master/graphics/lbm_taichi/lid_validation.png" height="293px">

## Example2: Kármán Vortex Street
<div align="center">
<img src="https://raw.githubusercontent.com/hietwll/common_files/master/graphics/lbm_taichi/VortexStreet.jpg" height="200px">
</div>

Kármán vortex street is an interesting phenomena in fluid dynamics. When fluids flow pass blunt body (say a cylinder), there exists a repeating pattern of swirling vortices, caused by a process known as vortex shedding. The Renolds number of this flow is defined as ``Re = U * D / niu``, where ``D`` means diameter. A solver with ``Re = 200`` can be given by:
```
lbm = lbm_solver(401, 101, 0.005, [0, 0, 1, 0],
      [[0.1, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
      1, [80.0, 50.0, 10.0])
```

<div align="center">
<img src="https://raw.githubusercontent.com/hietwll/common_files/master/graphics/lbm_taichi/karman.gif" height="150px">
</div>
