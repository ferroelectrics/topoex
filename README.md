# TOPO

Simulation of topological excitations in ferroic materials using [FEniCS](https://fenicsproject.org/) package.

## Usage notes

First of all you need to install FEniCS package. Please, follow [this](https://fenics.readthedocs.io/en/latest/installation.html) instructions.

Having the FEniCS installed and assuming that you have files 'meshgen.py', 'equations.py' and 'driver.py' in one directory with generated finite element mesh you can run the simulation using the following commands

```
python3 meshgen.py --msh_mesh your_mesh_name.msh --dim 3d # to generate mesh suitable for FEniCS
mpirun -n N_of_threads python3 driver.py --alpha1 -1.0 --alpha11 0.5 etc...
```
Use the values of Ginzburg-Landau expansion coefficients that suit the choosen material.

Please note that you need to generate the supporting mesh using for example the [gmsh](https://gmsh.info/) package.

Please, consult the documentation in tex file for greater level of technical details.
