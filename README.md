# 🧊 CrystalSAT  

**CrystalSAT** is a Python library for generating valid crystal structures using **constraint programming** and **SAT solvers**.  

Traditional crystal generators often rely on stochastic sampling or machine learning, which can produce large numbers of invalid candidates that require expensive filtering. CrystalSAT takes a different approach: by encoding chemical and geometric rules directly into a satisfiability framework, it produces **valid-by-construction** structures.  

CrystalSAT is designed for use in **materials science, solid-state physics, and computational chemistry**, providing a flexible and extensible foundation for research into lattice geometries, symmetry, and crystal structure prediction.  

---

## 🧰 Key Libraries  

| Library        | Purpose                                                                |
|----------------|------------------------------------------------------------------------|
| `pymatgen`     | Crystal manipulation (lattices, symmetry, CIF support)                 |
| `ASE`          | Atomic simulations and structure conversion/visualisation              |
| `pySAT`        | SAT/constraint solving for structure generation                        |
| `pypblib`      | Optimised backend for pseudo-Boolean constraints                       |
|  `NumPy`       | Library for mathematics calculations                                   |

---

## 📂 Project Structure  

```plaintext
crystalsat/
├── notebooks/       # Example notebooks and tutorials
├── experiments/     # Benchmarks and experimental runs
├── docs/            # Documentation and references
├── src/             # Core source code
│   └── constraints/ # Constraint encodings and solver logic
├── data/            # CIF files and generated structures
├── env/             # Conda environment setup
├── README.md        # Project overview
└── .gitignore       # Ignore cache/build files

