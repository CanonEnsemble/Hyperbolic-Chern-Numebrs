# Hyperbolic Chern Number
This repository contains code and data accompanying Ref. [1]. They are used to compute the Chern number for each irreducible representation (irrep) of the {8,3} Haldane model [2]. The project is divided into two main sections: Hamiltonian construction and Chern number calculation. 

## Construction of Hamiltonian 

The folder `0_ham` contains the codes and data necessary for generating the symbolic Hamiltonian for the {8,3} Haldane model. The symbolic Hamiltonian is generated using Python (version 3.6.9).

- Coset tables:

  - We have included the coset tables for five periodic clusters, `coset_table_G20_abelian.txt`, `coset_table_G24_nonabelian.txt`, `coset_table_G48_nonabelian.txt`, `coset_table_G56_nonabelian.txt`, and `coset_table_G100_nonabelian.txt`.
  - Each file is formatted as follows:
  The first two columns are indices for two sites $(i, j)$. The remaining eight columns contain binary values (1 or 0), indicating if $i$ and $j$ are related by one of the four generators (and their inverses).


- `0_cluster-Hh-flux.ipynb`:

  A Jupyter notebook that constructs the {8,3} Haldane model using coset tables.
  The output is a symbolic Hamiltonian saved as a pickle file.


- `1_pickle_to_mat.ipynb`:

  A Python script that converts the symbolic Hamiltonian from a pickle file to a MATLAB .mat file for further processing.
  The .mat file contains a four-component array `HAll`:
  The first two components represent the flux points $(\phi_i, \phi_j)$.
  The latter two components store the Hamiltonian matrix.


## Computation of Chern number
The files in the folder `1_chern` are used to compute the Chern number associated with a particular irrep. The codes are written in MATLAB, but the necessary input files are provided as .txt files if one wants to use a different language.

- Projector matrices:
  
  Inside the subfolders `G24` and `G48` are .txt files that give the projector matrices onto irreps of two translation groups of orders 24 and 48, respectively. Their character tables are listed here:
  
  
  For example, the files `projector_G24_1_real.txt` and `projector_G24_1_imag.txt` contain the real and imaginary parts, respectively, of the projector matrices associated with the irrep `1` in the group of order 24.

- `0_chern_number.m`:
  
  This code computes the Chern number associated with an irrep, taking the projector matrices and Hamiltonian as input. It uses the function `simdiag.m` to simultaneous diagonalize the projector and Hamiltonian [3]. 



## References
[1] C. Sun, A. Chen, T. Bzdušek, and J. Maciejko, *Topological linear response of hyperbolic Chern insulators*, [arXiv:2406.08388](https://doi.org/10.48550/arXiv.2406.08388), 2024

[2] D. Urwyler, P. Lenggenhager, I. Boettcher, R. Thomale, T. Neupert, and T. Bzdušek, *Hyperbolic Topological Band Insulators*,
[Phys. Rev. Lett. 129, 246402](https://doi.org/10.1103/PhysRevLett.129.246402), 2022

[3] Chris, *simdiag.m*, [MATLAB Central File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/46794-simdiag-m), version 1.0.0.0, 2024.

## Citation
If you use this code, please cite [this paper](https://doi.org/10.48550/arXiv.2406.08388):

```
@article{Sun2024Topological,
  title={Topological linear response of hyperbolic Chern insulators},
  author={Sun, Canon and Chen, Anffany and Bzdu{\v{s}}ek, Tom{\'a}{\v{s}} and Maciejko, Joseph},
  journal={arXiv preprint arXiv:2406.08388},
  year={2024},
  doi={https://doi.org/10.48550/arXiv.2406.08388}
}
