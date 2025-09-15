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
├── Base.py                 # Base classes and shared functionality
├── Cardinality.py          # Cardinality constraints (min/max atom counts, etc.)
├── Constraints.py          # Core constraint definitions
├── Coordinate.py           # Coordinate handling and transformations
├── Encoding.py             # CNF encodings and SAT solver interfaces
├── Get.py                  # Query helpers for retrieving constraints/data
├── Grab.py                 # Utility functions for input/output operations
├── NeighborAndDistances.py # Neighbor search and distance calculations
├── NeighborConstraints.py  # Constraints based on neighbor relations
├── OrbitsAndSymmetry.py    # Symmetry operations and orbit representations
├── SolveAndExport.py       # Running solvers and exporting valid structures
├── __init__.py             # Package initialisation
├── shannon-radii.json      # Ionic radii reference data
└── temp/                   # Temporary files or cached data


