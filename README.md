# TOPO

Simulation of topological excitations if ferroic materials using [FEniCS](https://fenicsproject.org/) package.

## Usage notes

First of all you need to install FEniCS package. Please, follow [this](https://fenics.readthedocs.io/en/latest/installation.html) instructions.

Having the FEniCS installed you can run the simulation using follow command

```
mpirun -n 8 python3 driver.py --alpha1 -1.0 --alpha11 0.5 etc...
```
