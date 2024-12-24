# Stokes Flow of an Evolving Fluid Film with Arbitrary Shape and Topology
This repository provides a minimal MATLAB implementation of a fluid membrane undergoing Helfrich-Stokes relaxation, as described in the paper:

**Stokes Flow of an Evolving Fluid Film with Arbitrary Shape and Topology**

*Cuncheng Zhu, David Saintillan, and Albert Chern*

## Getting Started
This implementation relies on the following key components:
* `Mesh.m`: Constructs the halfedge data structure for a triangular manifold mesh.
* `Geometry.m`: Computes geometric attributes of the $\mathbb{R}^3$-embedded mesh and the spatial differential operators for an evolving surface.
* `main.m`: The main executable that performs the variational time integration.

## Prerequisites
To run this code, you will need:
* MATLAB
* The [`sptensor`](https://www.tensortoolbox.org/sptensor_doc.html) library for efficient storage and manipulation of sparse multidimensional arrays.

## Running the Code
1. Clone this repository:
```Bash
git clone https://github.com/CunchengZhu/Evolving-Stokes-flow-2024.git
```
2.	Open `main.m` in MATLAB and execute it to start the simulation.

## Visualization
* **MATLAB Visualization**: Use the built-in `trisurf` function to visualize the triangular mesh and the `quiver` function to display the fluid velocity. 
* **External Tools**: For advanced visualization, export the mesh data and use tools like [Houdini](https://www.sidefx.com/download/). Learn more about Houdini with this [Introduction to Houdini](https://cseweb.ucsd.edu/~alchern/teaching/houdini/).

## Resources and Attribution
If you use this code in your academic projects or wish to learn more about the underlying methodology, please refer to and cite our papers:

### Main Paper
```bibtex
@article{zhu2024stokesflowevolvingfluid,
    title={Stokes flow of an evolving fluid film with arbitrary shape and topology},
    author={Zhu, Cuncheng and Saintillan, David and Chern, Albert},
    journal={arXiv preprint arXiv:2407.14025},
    year={2024}
}
```

### Helfrich Energy
The discretization of the Helfrich energy is based on the following paper: 

```bibtex
@article{Zhu_Mem3DG_Modeling_Membrane_2022,
    author = {Zhu, Cuncheng and Lee, Christopher T. and Rangamani, Padmini},
    doi = {10.1016/j.bpr.2022.100062},
    journal = {Biophysical Reports},
    title = {{Mem3DG: Modeling Membrane Mechanochemical Dynamics in 3D using Discrete Differential Geometry}},
    year = {2022}
}
```