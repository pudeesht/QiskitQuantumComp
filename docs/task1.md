# Task 1 Documentation: VQE for H₂ Molecule

## VQE Configuration

### Molecule Geometry, Basis Set, Charge, Spin
```python
molecule = MoleculeInfo(["H", "H"], [(0.0, 0.0, 0.0), (0.0, 0.0, r)], charge=0, multiplicity=1)
driver = PySCFDriver.from_molecule(molecule, basis="sto3g")
problem = driver.run()
```
- **Molecule**: Hydrogen (H₂)
- **Geometry**: Two hydrogen atoms with varying bond distances between 0.2 Å and 2.0 Å.
- **Basis Set**: STO-3G
- **Charge**: 0
- **Multiplicity (Spin)**: 1 (Singlet state)

### Second Quantized Hamiltonian
The molecular electronic structure problem is encoded using Qiskit Nature’s `PySCFDriver`, which internally uses PySCF to compute the integrals and returns a `ElectronicStructureProblem`. The Hamiltonian is constructed in its second quantized form.


### Mapper
```python
mapper = JordanWignerMapper()
```
- **Mapper Used**: Jordan-Wigner Mapper
- This maps the fermionic second quantized operators to qubit operators compatible with quantum circuits.

### Ansatz
```python
ansatz = UCCSD(
    problem.num_spatial_orbitals,
    problem.num_particles,
    mapper,
    initial_state=HartreeFock(
        problem.num_spatial_orbitals,
        problem.num_particles,
        mapper
    )
)
```
- **Ansatz Used**: Unitary Coupled Cluster Singles and Doubles (UCCSD)
- **Initial State**: Hartree-Fock
- The Hartree-Fock state is constructed using the number of spatial orbitals and particles derived from the molecule.

### Optimizer
- **Optimizer**: SLSQP (Sequential Least Squares Programming)
- A gradient-based classical optimizer used to minimize the expectation value of the energy.
- Configured with `maxiter=1000` for better convergence.

### Backend / Estimator
- **Estimator Used**: AerEstimator from Qiskit Aer (simulator)
- Used with 1024 shots and a fixed random seed (42) to ensure reproducibility.

### Initial Parameters
- The initial parameters for the ansatz circuit are set to zeros using:
  ```python
  vqe.initial_point = np.zeros(ansatz.num_parameters)
  ```

## Description of Results Obtained

This script performs a **Potential Energy Surface (PES)** scan for the H₂ molecule by varying the bond length from 0.2 Å to 2.0 Å in 20 steps. For each geometry, it:
1. Constructs the second quantized Hamiltonian using PySCF.
2. Maps it to a qubit operator using the Jordan-Wigner transform.
3. Prepares a UCCSD ansatz with a Hartree-Fock initial state.
4. Uses the SLSQP optimizer with an AerEstimator backend to perform VQE and estimate the ground state energy.

The results are plotted using matplotlib to visualize the relationship between bond length and ground state energy, giving a curve characteristic of the H₂ potential energy surface.