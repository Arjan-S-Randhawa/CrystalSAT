# 🧊 Crystal Learning

This repository documents my learning journey and experimentation with materials science, structure prediction, and symbolic AI tools — particularly for my Google DeepMind internship. It focuses on using Python libraries such as `pymatgen`, `ASE`, and `python-sat` to explore, simulate, and reason about atomic structures.

---

## 🧪 Goals

- Learn to manipulate and analyse **crystal structures**
- Simulate lattice geometries and calculate material properties
- Explore **constraint-solving and symbolic reasoning** using SAT solvers
- Build clean and reproducible environments with **Conda**
- Document insights and methods for future academic work

---

## 🧰 Key Libraries Used

| Library        | Purpose                                                                |
|----------------|------------------------------------------------------------------------|
| `pymatgen`     | Materials structure manipulation (lattices, unit cells, symmetry)      |
| `ASE`          | Atomic simulations (building, viewing, converting structures)          |
| `python-sat`   | Constraint solving using SAT solvers (e.g., for symbolic reasoning)    |
| `pypblib`      | Optional backend for `python-sat` (optimised for pseudo-Boolean logic) |

---

## 📂 Project Structure
crystal-learning/
├── notebooks/         # Jupyter notebooks with live experiments
├── experiments/       # Writeups and notes on experiment results
├── explanations/      # Markdown explanations of source code + key formulas
├── src/               # Main Python source code
│   └── utils/         # Helper functions for I/O, visualisation, etc.
├── data/              # Test CIFs or simulation files
├── env/               # environment.yml file (Conda environment setup)
├── README.md          # Project overview and documentation
└── .gitignore         # Prevents pushing unnecessary files to GitHub
