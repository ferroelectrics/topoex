# TOPO

Simulation of topological excitations if ferroic materials using [FEniCS](https://fenicsproject.org/) package.

## Usage notes

First of all you need to install FEniCS package. Please, follow [this](https://fenics.readthedocs.io/en/latest/installation.html) instructions.

Having the FEniCS installed you can run the simulation using follow command

```
mpirun -n N_of_threads python3 driver.py --alpha1 -1.0 --alpha11 0.5 etc...
```
Use the values of Ginzburg-Landau expansion that suit the choosen material.

Please note that you need to generate the supporting mesh using for example the [gmsh](https://gmsh.info/) package.
