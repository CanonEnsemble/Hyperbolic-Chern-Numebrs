#-----------------------------------------------------------------------------------
#   cluster-Hh-flux.py
# 
#   Construct the Haldane model (with flux threading) on a {8,3} PBC cluster built 
#   from the coset table of a given {8,8} PBC cluster. 
#
#   Code written by Anffany Chen. Last updated 2024-03-18.
#-----------------------------------------------------------------------------------
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import sympy as sp
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct
import time
import pickle

# Path to the coset table of a given {8,8} PBC cluster
infile = './data/coset_table_G24_nonabelian.txt'

# Path for saving the resulting Haldane-flux Hamiltonian
outfile = './data/Hh-flux-G24.pickle'

# Haldane model parameters
t1 = 1
t2 = 1/6
mass = 1/3

def main(): 
    start = time.time()

    # Construct hopping matrices
    T1, V1 = nn_hoppings()
    T2, T3, V2 = nnn_hoppings(f = 1j)

    # Construct matrix for sublattice mass 
    M = sublattice_mass()

    # Import {8,8} cluster 
    amat = LoadAdjacencyMatrix(infile)
    nn_gens = LoadGeneratorData(infile)

    # Define symbols for flux phases 
    phi = sp.symbols(f"phi0:{8}", real=True)    

    # Build {8,3} Haldane model (with flux threading) on the {8,8} PBC cluster 
    total = len(amat[0,:])
    H1 = TensorProduct(sp.Matrix(np.eye(total)), sp.Matrix(V1))
    H2 = TensorProduct(sp.Matrix(np.eye(total)), sp.Matrix(V2))
    for n in range(total):
        for m in range(total):
            if amat[n,m] > 1:
                print('Error: This cluster has a multi-bond!')
            elif amat[n,m] == 1:
                tempM = sp.zeros(total,total)
                tempM[n,m] = 1
                gen_idx = np.nonzero(np.squeeze(nn_gens[n,m,:]))[0][0]   
                H1 += TensorProduct(tempM, sp.exp(1j*phi[gen_idx])*sp.Matrix(np.squeeze(T1[:,:,gen_idx])))
                H2 += TensorProduct(tempM, sp.exp(1j*phi[gen_idx])*sp.Matrix(np.squeeze(T2[:,:,gen_idx])))

                # Implement nnn hoppings that are two unit cells apart
                m_NN = np.nonzero(amat[:,m])[0]
                for m_nn in m_NN: 
                    second_gen_idx = np.nonzero(np.squeeze(nn_gens[m,m_nn,:]))[0][0]
                    if second_gen_idx == (gen_idx+3)%8: 
                        tempM = sp.zeros(total,total)
                        tempM[n,m_nn] = 1
                        H2 += TensorProduct(tempM, sp.exp(1j*(phi[gen_idx]+phi[second_gen_idx]))*sp.Matrix(np.squeeze(T3[:,:,gen_idx])))
                    elif second_gen_idx == (gen_idx+5)%8: 
                        tempM = sp.zeros(total,total)
                        tempM[n,m_nn] = 1
                        H2 += TensorProduct(tempM, sp.exp(1j*(phi[gen_idx]+phi[second_gen_idx]))*sp.Matrix(np.squeeze(T3[:,:,(gen_idx+1)%8]).T.conj()))

    # Reduce to k1, k2, k3, k4
    H1 = H1.subs([(phi[4],-phi[0]),(phi[5],-phi[1]),(phi[6],-phi[2]),(phi[7],-phi[3])])
    H2 = H2.subs([(phi[4],-phi[0]),(phi[5],-phi[1]),(phi[6],-phi[2]),(phi[7],-phi[3])])

    # Check Hermiticity    
    print('Checking hermiticity of H2...')
    print(np.unique(np.abs(np.array(H2 - Dagger(H2)).astype(complex)), return_counts=True))

    # Check connectivity
    H1_noflux = np.array(H1.subs([(phi[0],0),(phi[1],0),(phi[2],0),(phi[3],0)])).astype(float)
    H2_noflux = np.array(H2.subs([(phi[0],0),(phi[1],0),(phi[2],0),(phi[3],0)])).astype(complex)
    print('Checking connectivity of H2 (without flux threading)...')
    print(np.unique(H1_noflux@H1_noflux - 3*np.eye(total*16) - np.abs(H2_noflux), return_counts=True))

    H = t1*H1 + t2*H2 + mass*TensorProduct(sp.Matrix(np.eye(total)), sp.Matrix(M))

    # Pickle to file
    with open(outfile, 'wb') as outf:
        outf.write(pickle.dumps(H))

    # Plot the eigenvalues at zero flux 
    eval = LA.eigvalsh(np.array(H.subs([(phi[0],0),(phi[1],0),(phi[2],0),(phi[3],0)])).astype(np.complex128))
    plt.plot(np.arange(1, total*16+1), eval, 'g.', ms = 0.5)
    plt.plot(np.arange(1, total*5+1), eval[:int(total*5)], 'r.', ms = 0.5)
    plt.plot(np.arange(total*11+1,total*16+1), eval[-int(total*5):], 'r.', ms = 0.5)
    plt.xlabel('Sorted Eigenvector Index')
    plt.ylabel('Eigenvalue')
    plt.title(f'{{8,3}} Haldane model with {total:.0f} unit cells')
    plt.show()

    print(f'Total time lapsed: {time.time()-start:.0f} seconds.')

def nn_hoppings():
    T1 = np.zeros((16,16,8))
    V1 = np.zeros((16,16))
    for a in range(8):  
        T1[8+a, 8+(a+5)%8, a] = 1
        if a < 4: 
            T1[8+(a+5)%8, 8+a, a+4] = 1
        else: 
            T1[8+(a+5)%8, 8+a, a%4] = 1
        V1[(a+1)%8, a] = 1
        V1[a+8, a] = 1
    V1 = V1 + V1.T
    return T1, V1

def nnn_hoppings(f):
    T2 = np.zeros((16,16,8), dtype=complex)  # nnn hoppings between the central polygon and a neighboring polygon
    T3 = np.zeros((16,16,8), dtype=complex)  # nnn hoppings that are two unit cells apart
    V2 = np.zeros((16,16), dtype=complex)    # nnn hoppings inside the central polygon

    for a in range(8):
        V2[(a+2)%8, a] = f
        V2[a, 8+(a+1)%8] = f
        V2[(a-1)%8+8, a] = f

        T2[a, (a+5)%8+8, a] = -f
        if a < 4:
            T2[a+1, a+12, a] = f
        else:
            T2[(a+1)%8, a+4, a] = f

        T2[a+8, (a+5)%8, a] = -f
        if a < 7: 
            T2[a+9, (a+4)%8, a] = f
        else:
            T2[a+1, (a+4)%8, a] = f

        T3[(a+3)%8+8, (a+5)%8+8, a] = f  
        
    V2 = V2 + V2.T.conj()
    return T2, T3, V2

def sublattice_mass():
    M = np.zeros((16, 16))
    for a in range(8): 
        M[a,a] = (-1)**a
        M[a+8, a+8] = -(-1)**a    
    return M

def LoadAdjacencyMatrix(filename):
    data= np.fromfile(filename, dtype=int, sep=" ").reshape(-1,10)
    row = data[:,0] - 1 # python indices begin with 0
    col = data[:,1] - 1
    gens = data[:,2:10]
    sA = coo_matrix((np.sum(gens,axis=1),(row,col))).tocsr()
    return np.asarray(sA.todense())

def LoadGeneratorData(filename):
    data= np.fromfile(filename, dtype=int, sep=" ").reshape(-1,10)
    row = data[:,0] - 1 # python indices begin with 0
    col = data[:,1] - 1
    gens = data[:,2:10]
    total = max(row) + 1
    genMatrix = np.zeros((total, total , 8))
    for i in range(len(row)):
        site1 = row[i]
        site2 = col[i]
        genMatrix[site1, site2, :] = gens[i,:]
    return genMatrix

if __name__ == '__main__':
    main()
