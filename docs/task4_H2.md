# 🧪 Comparison of H₂ Ground State Energy Calculations

This document summarizes the computed ground state energies for the H₂ molecule using three different quantum chemistry approaches: classical CCSD via PySCF, VQE without circuit optimizations, and VQE with FreezeCore and qubit tapering.

---

## ⚙️ Common Setup
- **Molecular Geometry**: H – H bond distance = 0.735 Å  
- **Basis Set**: STO-3G  
- **Charge**: 0  
- **Multiplicity**: 1  

---

## 1️⃣ PySCF: CCSD (Classical Reference)

| Method | Description | Energy (Hartree) |
|--------|-------------|------------------|
| PySCF | RHF → CCSD | `mycc.e_tot` |

- Used PySCF to perform Restricted Hartree-Fock (RHF) followed by CCSD.  
- Acts as the classical benchmark for accuracy.

---

## 2️⃣ VQE (No FreezeCore, No Tapering)

| Component | Value |
|----------|--------|
| Ansatz | UCCSD |
| Mapper | Jordan-Wigner |
| Optimizer | COBYLA |
| Initial Parameters | All Zeros |
| Optimization Level | 1 (transpiler) |

- **Nuclear repulsion energy**: `problem.nuclear_repulsion_energy`  
- **Final VQE result**: `result.fun`  
- **Final Total Energy**: `result.fun + nuclear_repulsion_energy`  

This is a full-space simulation without any optimization in terms of qubit reduction.

---

## 3️⃣ VQE (With FreezeCore + Qubit Tapering)

| Component | Value |
|----------|--------|
| Ansatz | UCCSD |
| Mapper | Jordan-Wigner with qubit tapering |
| Optimizer | COBYLA |
| Initial Parameters | All Zeros |
| Optimization Level | 1 (transpiler) |

- **Frozen core applied** using `FreezeCoreTransformer`  
- **Tapering applied** via `get_tapered_mapper(...)`  
- **Frozen Core Energy Constant**: `transformed_problem.hamiltonian.constants.get("FreezeCoreTransformer", 0.0)`  
- **Final VQE result**: `result.fun`  
- **Final Total Energy**: `result.fun + nuclear_repulsion_energy + frozen_core_energy`  

---



## 📝 Notes

- All simulations use **EstimatorV2** and **UCCSD** as the ansatz for fair comparison.  
- Optimization history (energy per iteration) can also be visualized using `cost_history`.

---

## Result Comparison Summary

| Method                                | Total Energy (Hartree)     | Notes                                                   |
|---------------------------------------|-----------------------------|---------------------------------------------------------|
| **H₂ - PySCF CCSD**                   | -1.1373061933919688         | High-accuracy classical benchmark using PySCF           |
| **H₂ - VQE (no tapering)**            | -1.137306                   | GroundState: -1.857275 + Nuclear: 0.719969              |
| **H₂ - VQE (tapering + freeze core)** | -1.137306                   | GroundState: -1.857275 + Nuclear: 0.719969 + FC: 0.0    |
