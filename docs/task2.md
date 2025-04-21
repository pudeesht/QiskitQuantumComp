
# Task 2: VQE using EstimatorV2, Custom Loop, and Transpilation Analysis

## Objective

This task demonstrates a lower-level manual implementation of the VQE algorithm using `EstimatorV2`, bypassing the built-in `VQE` class. We also analyze how transpilation optimization levels affect circuit complexity.

### Goals
1. Use `EstimatorV2` from `qiskit_aer.primitives` to compute expectation values.
2. Run VQE manually by defining a cost function and using a classical optimizer.
3. Transpile the ansatz circuit at optimization levels 0–3 using `preset_pass_manager` and benchmark the circuit's depth and width.
4. Track the convergence of the energy for each transpilation level.

## Configuration and Execution

### Molecule and Hamiltonian

First, we define the molecular system. Here, we use the H₂ molecule with a bond length of 0.735 Å and the STO-3G basis set.

```python
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper

driver = PySCFDriver(atom="H 0 0 0; H 0 0 0.735", basis="sto3g")
problem = driver.run()
hamiltonian = problem.hamiltonian.second_q_op()
mapper = JordanWignerMapper()
qubit_op = mapper.map(hamiltonian)
```

- **PySCFDriver**: The `PySCFDriver` from `qiskit_nature` is used to generate the molecular Hamiltonian. It takes in the atomic structure and basis set to create a quantum chemistry problem.
- **JordanWignerMapper**: This is the mapping used to convert fermionic operators (from the molecular Hamiltonian) into qubit operators, which is necessary for quantum computation.

### Ansatz and Initial State

The next step is to define the ansatz (variational form) for the VQE. We use **UCCSD** (Unitary Coupled Cluster Singles and Doubles) as the ansatz, with the Hartree-Fock state as the initial state.

```python
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD

num_particles = problem.num_particles
num_spatial_orbitals = problem.num_spatial_orbitals
hf_state = HartreeFock(num_spatial_orbitals, num_particles, mapper)
ansatz = UCCSD(num_spatial_orbitals, num_particles, mapper, initial_state=hf_state)
```

- **HartreeFock**: The `HartreeFock` method is used to create the initial state for the quantum system before applying the UCCSD operator. 
- **UCCSD**: The Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz is widely used in quantum chemistry for variational quantum eigensolvers (VQE), capturing the correlation between electrons in a molecule.

### Backend and Estimator

We use the `BasicSimulator` from Qiskit’s Aer provider as the backend for simulation, and the `EstimatorV2` from `qiskit_aer.primitives` to evaluate expectation values.

```python
from qiskit_aer.primitives import EstimatorV2
from qiskit.providers.basic_provider import BasicSimulator

backend = BasicSimulator()
estimator = EstimatorV2()
```

- **EstimatorV2**: This is a quantum primitive used to run quantum circuits and measure expectation values. It provides a more efficient and flexible way of evaluating quantum circuits compared to the older `Estimator` class.

## Transpilation and Manual VQE Loop

The core of the task involves transpiling the ansatz circuit at different optimization levels, running VQE manually, and tracking the convergence.

We iterate over the four optimization levels: 0 (no optimization), 1, 2, and 3 (maximum optimization). At each level, we transpile the ansatz circuit using Qiskit’s `preset_pass_manager` and benchmark its depth and width. The **depth** of a circuit refers to the longest path from the input to the output (number of gates in the longest sequence), and the **width** refers to the number of qubits used.

```python
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from scipy.optimize import minimize
import numpy as np

convergence_data = {}

for level in range(4):
    print(f"
 Optimization Level {level} ")

    pm = generate_preset_pass_manager(backend=backend, optimization_level=level)
    transpiled_ansatz = pm.run(ansatz)
    initial_params = np.zeros(transpiled_ansatz.num_parameters)

    print(f"Initial params {initial_params} ")
    print("Circuit depth:", transpiled_ansatz.depth())
    print("Circuit width:", transpiled_ansatz.width())

    observable = qubit_op.apply_layout(transpiled_ansatz.layout)
    cost_history_dict = {"iters": 0, "cost_history": []}

    def cost_func(params):
        pub = (transpiled_ansatz, [observable], [params])
        result = estimator.run([pub]).result()
        energy = result[0].data.evs[0]
        cost_history_dict["iters"] += 1
        cost_history_dict["cost_history"].append(energy)
        print(f"  Iteration {cost_history_dict['iters']} -> Energy: {energy:.6f}")
        return energy

    result = minimize(cost_func, initial_params, method="COBYLA")
    convergence_data[level] = cost_history_dict["cost_history"]
```

### Explanation of the Code:

1. **Preset Pass Manager**: 
   The `generate_preset_pass_manager` function is used to generate a pass manager with different optimization levels (0-3). These levels define how much optimization should be performed on the quantum circuit during transpilation. Higher levels result in more aggressive optimizations, aiming to reduce circuit depth and width.
   
2. **Transpiled Circuit**:
   The ansatz circuit is transpiled according to the optimization level, and the depth and width of the transpiled circuit are printed. These values give us an idea of how the transpiler affects the circuit complexity at each optimization level.

3. **Cost Function**:
   A cost function is defined for the VQE. It takes the current parameter values, applies them to the ansatz, and computes the expectation value (energy) using the `EstimatorV2` object.

4. **Optimizer**:
   The optimizer used in this case is COBYLA (a non-gradient-based method). The `minimize` function from `scipy.optimize` is used to perform the optimization loop, minimizing the cost function.

5. **Convergence Tracking**:
   We track the convergence of the optimization by storing the energy values in `cost_history_dict["cost_history"]` for each iteration.

## Analysis and Observations

- **Transpilation Impact**:
  - As the optimization level increases, the depth and width of the quantum circuits generally decrease. This shows that the transpiler is able to optimize the circuit by reducing redundant operations and simplifying gate sequences.
  - The higher optimization levels may also result in fewer gates and thus, potentially faster execution on real quantum hardware.

- **Manual VQE Execution**:
  - Instead of using the built-in `VQE` class from `qiskit_algorithms`, we manually define a cost function to compute expected energy values for each iteration. This gives us more flexibility in defining the optimization process.

- **Tracking Convergence**:
  - The convergence of the energy is tracked for each optimization level and stored in `convergence_data`. This allows us to compare the performance of the VQE optimization process at different transpilation levels.

- **Circuit Complexity**:
  - By checking the circuit depth and width at different optimization levels, we can analyze how efficient the transpiler is at optimizing the quantum circuit.

---

This task provides insight into the lower-level workings of the VQE algorithm and how transpilation optimization levels can affect circuit complexity and performance. It also shows how to manually perform VQE without relying on the `VQE` class, offering more control over the optimization process.
