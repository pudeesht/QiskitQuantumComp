
#  Task 3: Freeze Core & Qubit Tapering on LiH Molecule

##  Objective  
To simplify the quantum simulation of the **LiH molecule** using two techniques:  
- **Freeze Core Approximation**: Removes core (inert) electrons from the active space.  
- **Qubit Tapering**: Uses symmetries in the Hamiltonian to reduce qubit count.  

The aim is to preserve chemical accuracy while reducing simulation complexity.

---

##  Code Explanation

### 🔹 Molecule & Problem Setup

```python
molecule = MoleculeInfo(["Li", "H"], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.6)], charge=0, multiplicity=1)
driver = PySCFDriver.from_molecule(molecule, basis="sto3g")
problem = driver.run()
```

- **MoleculeInfo** defines LiH in 3D space.
- `PySCFDriver` computes the electronic structure.
- `problem` is the full electronic structure problem (Hamiltonian, orbitals, particles, etc.).

---

### 🔹 Applying Freeze Core Approximation

```python
fc_transformer = FreezeCoreTransformer()
transformed_problem = fc_transformer.transform(problem)
```

- **What it does**: Removes core orbitals (1s² of Li) that don’t participate in bonding.  
- **Why**: These core electrons don’t affect the chemical properties much, so excluding them reduces problem size.

 **Effect**:
- Orbitals drop from 12 → 10
- Electrons drop from (2, 2) → (1, 1)

---

### 🔹 Qubit Mapping via Jordan-Wigner

```python
qubit_op = JordanWignerMapper().map(transformed_problem.second_q_ops()[0])
hamiltonian = SparsePauliOp(qubit_op)
```

- Converts the fermionic Hamiltonian into qubit operators using **Jordan-Wigner mapping**, which translates creation/annihilation operators to Pauli strings.
- `SparsePauliOp` wraps this in a more efficient format for manipulation.

---

### 🔹 Finding Z₂ Symmetries & Tapering Qubits

```python
z2_symmetries = Z2Symmetries.find_z2_symmetries(hamiltonian)
```

- This scans the Hamiltonian for **Z₂ symmetries** — Pauli operators that commute with the full Hamiltonian and square to identity.
- Each such symmetry lets us **eliminate one qubit** via **tapering**.

If symmetries are found, we taper:

```python
tapering_values = [1 if i % 2 == 0 else -1 for i in range(len(z2_symmetries.symmetries))]
z2_symmetries.tapering_values = tapering_values
tapered_hamiltonian = z2_symmetries.taper(hamiltonian)
```

- Tapering values manually assigned (could be chosen based on known spin states).
- Applies the tapering, producing a **simpler Hamiltonian** with fewer qubits.

 **Effect**:  
- Qubit count reduces from 10 → **6**

---

### 🔹 Output & Interpretation

```python
print("Before FreezeCore:", problem.num_spin_orbitals, problem.num_particles)
print("After FreezeCore:", transformed_problem.num_spin_orbitals, transformed_problem.num_particles)
print("Original Hamiltonian:", hamiltonian)
print("Tapered Hamiltonian:", tapered_hamiltonian)
```

This verifies:
- FreezeCore reduced orbitals and electrons
- Tapering reduced the qubit Hamiltonian complexity
- Final Hamiltonian is significantly easier to simulate on quantum hardware

---

##  Summary Table

| Stage              | Orbitals | Particles | Qubits |
|--------------------|----------|-----------|--------|
| Original (LiH)     | 12       | (2, 2)    | 12     |
| After Freeze Core  | 10       | (1, 1)    | 10     |
| After Tapering     | —        | —         | **6**  |

---

##  Conclusion

This task demonstrates how **quantum chemistry preprocessing** techniques like FreezeCore and Qubit Tapering can:
- Simplify Hamiltonians,
- Reduce resource requirements,
- Make quantum simulations more tractable.

Both steps preserve core physics while cutting down on qubits — a vital step for running quantum algorithms on near-term devices.
