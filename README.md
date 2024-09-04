# Hyperbolic Chern Number
This repository contains code and data for computing the Chern number of the {8,3} Haldane model for each irreducible representation (irrep). The project is divided into two main sections: Hamiltonian construction and Chern number calculation. 

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
