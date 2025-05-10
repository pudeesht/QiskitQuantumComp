
# LiH - Quantum Chemistry Calculations with PySCF and VQE

## Overview
This document compares three different approaches to calculating the ground-state energy of the LiH molecule using quantum chemistry methods: 

1. **LiH - PySCF CCSD** (Standard PySCF method without FreezeCore)
2. **LiH - PySCF CCSD with FreezeCore** (Using the FreezeCore approximation)
3. **LiH - VQE with Tapering and FreezeCore** (Using the Variational Quantum Eigensolver (VQE) with Qiskit)

The LiH molecule has been optimized in all three cases, and the energy results have been compared. Below are the details for each case.

## 1. LiH - PySCF CCSD (No FreezeCore)

```python
from pyscf import gto, scf, cc

# Defining the LiH molecule
mol = gto.M(
    atom='Li 0 0 0; H 0 0 1.6',
    basis='sto-3g',
    unit='Angstrom',
    charge=0,
    spin=0,
)

# Performing RHF SCF calculation
mf = scf.RHF(mol).run()

# Running CCSD calculation
mycc = cc.CCSD(mf).run()

# Output CCSD energy
print("LiH PySCF CCSD Energy (no freeze core):", mycc.e_tot)
```

- **Result**: `LiH PySCF(CCSD) Energy: -7.882313821664073`

In this case, we use the coupled-cluster singles and doubles (CCSD) method in PySCF to calculate the energy of LiH without any approximations. The method computes the ground-state energy by iteratively solving for the wavefunction.

## 2. LiH - PySCF CCSD with FreezeCore

```python
mol = gto.M(
    atom='Li 0 0 0; H 0 0 1.6',
    basis='sto-3g',
    unit='Angstrom',
    charge=0,
    spin=0,
)

# Perform RHF SCF calculation
mf = scf.RHF(mol).run()

# Run CCSD calculation with FreezeCore
mycc = cc.CCSD(mf, frozen=1).run()

# Output CCSD energy with freeze core approximation
print("LiH PySCF CCSD Energy (with freeze core):", mycc.e_tot)
```

- **Result**: `LiH PySCF(CCSD) + FreezeCore Energy: -7.882096600731343`

In this case, the FreezeCore approximation is applied. FreezeCore essentially freezes the core orbitals (those that are far from the Fermi level) and treats them as if they do not interact in the correlation calculation. This reduces the size of the active space, speeding up the computation.

## 3. LiH - VQE with Tapering and FreezeCore

```python
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.formats.molecule_info import MoleculeInfo
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
from qiskit_aer.primitives import EstimatorV2
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit.providers.basic_provider import BasicSimulator
from scipy.optimize import minimize
import numpy as np

# Defining the molecule
molecule = MoleculeInfo(["Li", "H"], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.6)], charge=0, multiplicity=1)
driver = PySCFDriver.from_molecule(molecule, basis="sto3g")
problem = driver.run()

# Apply FreezeCore approximation
fc_transformer = FreezeCoreTransformer()
transformed_problem = fc_transformer.transform(problem)

# Obtain the Hamiltonian after FreezeCore
hamiltonian = transformed_problem.second_q_ops()[0]

# Map the Hamiltonian to qubits using Jordan-Wigner mapping
mapper = JordanWignerMapper()
tapered_mapper = transformed_problem.get_tapered_mapper(mapper)
qubit_op = tapered_mapper.map(hamiltonian)

# Setup backend and ansatz
backend = BasicSimulator()
num_particles = transformed_problem.num_particles
num_spatial_orbitals = transformed_problem.num_spatial_orbitals
hf_state = HartreeFock(num_spatial_orbitals, num_particles, tapered_mapper)
ansatz = UCCSD(num_spatial_orbitals, num_particles, tapered_mapper, initial_state=hf_state)

# Set initial parameters and optimizer
initial_params = np.random.uniform(-0.1, 0.1, ansatz.num_parameters)
estimator = EstimatorV2()
cost_history = []

def cost_func(params):
    bound_circ = ansatz.assign_parameters(params)
    pm = generate_preset_pass_manager(backend=backend, optimization_level=1)
    transpiled_ansatz = pm.run(bound_circ)
    pub = (transpiled_ansatz, [qubit_op], [np.array([])])
    result = estimator.run([pub]).result()
    energy = result[0].data.evs[0]
    cost_history.append(energy)
    return energy

# Run VQE optimization
result = minimize(cost_func, initial_params, method="COBYLA", options={'maxiter': 200})

# Print final results
vqe_energy = result.fun
nuclearenergy = problem.nuclear_repulsion_energy
frozen_core_energy = transformed_problem.hamiltonian.constants.get("FreezeCoreTransformer", 0.0)
total_energy = vqe_energy + frozen_core_energy + nuclearenergy

print(f"VQE energy (active space only): {vqe_energy:.6f} Hartree")
print(f"Nuclear repulsion energy: {nuclearenergy} Hartree")
print(f"Frozen core energy from hamiltonian constants: {frozen_core_energy} Hartree")
print(f"Total ground state energy of LiH: {total_energy:.6f} Hartree")
```

- **Result**: `LiH - with Tapering and FrozenCore (VQE Sim): -7.882097`
    - **Ground State Energy**: -1.078084
    - **FreezeCore Energy**: -7.796219568777053
    - **Nuclear Repulsion Energy**: 0.992207270475

In this case, we use the VQE with both FreezeCore and qubit tapering techniques. The FreezeCore approximation is applied in the same way as the previous case, and qubit tapering is used to reduce the number of qubits needed for the calculation. The VQE method finds the best parameters for the ansatz (UCCSD) that minimize the energy.

## Results Comparison

| **Method**                             | **Energy (Hartree)**          |**Notes**                                                                         |
|----------------------------------------|------------------------------|------------------------------------------------------------------------------------|
| **LiH - PySCF (CCSD)**                 | -7.882313821664073           |High-accuracy classical benchmark using PySCF                                       |
| **LiH - PySCF (CCSD + FreezeCore)**    | -7.882096600731343           |High-accuracy classical benchmark using PySCF                                       |
| **LiH - VQE with Tapering and FreezeCore** | -7.882097                |(GroundState: -1.078084 + Freezecore: -7.796219568777053 + Nuclear:0.992207270475)  |

### Final Remarks

The VQE method with FreezeCore and qubit tapering closely matches the results from the CCSD calculation (with and without FreezeCore) for the LiH molecule. The computational efficiency gained through these approximations (FreezeCore and qubit tapering) allows for a faster calculation without a significant loss of accuracy.

