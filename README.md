# Physics-Simulations
Here I have collected some of the numerous research side projects I did throughout my career as a curious high school/university student.

Feel free to contact me about these projects at misa.toman97@gmail.com and I will happily answer!

## Solar wind simulation.

The goal of this simulation is to convince me that the typical illustration of the Aurora Borealis explanation is correct. What is typically drawn is a particle revolving around the magnetic field lines in a helical motion. To demonstrate this, an initial volley of charged particles is launched at a magnetic dipole and the trajectory of these particles is plotted. We see that typical helical motion around magnetic field lines is seen but is highly squished.

## 1D Time independent Schrodinger equation solver.
For any kind of (well behaved, bounded) potential V(x). An iterative method is used to diagonalize the matrix of the discretized hamiltonian. Eigenvalues and eigenvectors are plotted along with the potential. Besides that, I also calculate the error in the estimate of the energy levels.

## Symmetry adapted basis of a group representation D(G) on vector space V.
Mathematica code that given the representation of a point group D(G) on some vector space V outputs the symmetry adapted basis for this vector space. A few applications of this are given.
### Molecular orbitals
Because the hamiltonian H commutes with the representation D(G) they have shared eigenspaces. This can be used to simplify the eigenproblem of the hamiltonian because Hamiltonian can be reduced to each irreducible space of the representation. To find the irreducible spaces the Group Projection Operator technique is employed. I give two concrete examples.
#### Water molecule molecular orbitals
are constructed from atomic orbitals, concretely two S orbitals on each hydrogen atom and three P orbitals on oxygen. The hamiltonian is reduced in the irreducible eigenspaces of the representation D(G) and its spectrum is calculated.
#### Benzene molecule molecular orbitals
are constructed from atomic orbitals, concretely six P_z orbitals on each carbon atom. The hamiltonian is reduced in the irreducible eigenspaces of the representation D(G) and its spectrum is calculated. For each example I also plot what these orbitals look like!
### Normal modes
The second example involves the vibrations of a molecule. The linearized equations of motion for a system of particles that moves around the equilibrium configuration can be written as x'' = -K.x, where x represents a vector made of degrees of freedom for each particle, '' is the second time derivative and K is the so called "dynamic matrix". The dynamic matrix is encoding the mutual Hook law-esque couplings between the particles. Upon making the anzac for normal modes, x = exp(-iwt) c, this turns the problem into an eigenproblem w^2 c = K.c, therefore, the frequency of oscillations are precisely the square roots of the eigenvalues of K while c are the eigenvectors of K.
If the system possesses a symmetry G then one can construct the so called Dynamic representation. It is made as a tensor product of the defining representation of SO(3) and the permutation representation over the atoms. The group elements rotate/reflect the oscillation directions and permute the atoms of the system. This representation commutes with K and its symmetry adapted basis represents the normal modes.
