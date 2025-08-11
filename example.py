# example_basic.py
# ----------------
# Minimal CrystalSAT usage with solve_multiple.
from CrystalSAT import *

# 1) Grid & cell parameters
n_x, n_y, n_z = 7, 7, 7                # grid sites
a, b, c = 17.0, 17.0, 17.0              # cell lengths (Å)
alpha, beta, gamma = 90.0, 90.0, 90.0   # angles (deg)

# 2) Allowed species (atoms or (symbol, charge, CN))
#    if "Allowed is empty, then ALL atoms + ions are used
allowed = [
    "Li", "Fe", "C", "Ti", "O", "H",
    ("Ca", 2, "VIII"), # ion
    ("Ru", 8, "IV"),   # ion
]

# 3) Create model
crystal = CrystalSAT(n_x, n_y, n_z, a, b, c, alpha, beta, gamma, allowed)

# 4) Basic constraints
# This should always be called first to set up the CNF and variable pool.
crystal.initialise(pack=True)   # one per site + no overlaps


#---- Helper functions ----#
# 1) Get atom ID
# (returns integer ID for atom or ion)
crystal.atom_id("Li")                # returns ID for Li atom
crystal.atom_id("Ca", 2, "VIII") # returns ID for Ca2+ ion (VIII CN)

# 2) Get inverse ID
# (returns symbol or (symbol, charge, CN) for atom or ion)
crystal.inverse_id(0)               # returns "Li" for Li atom
crystal.inverse_id(6)               # returns ("Ca", 2, "VIII") for Ca2+ ion (VIII CN)

# 3) manually encoding a variable:
# this gets the cnf variable ID for Li at position (0,0,0)
li_000 = crystal.encode_var(0,0,0, crystal.atom_id("Li"))
# 4) Append to CNF
# this allows user to create manual custom constraints
crystal.cnf.append(li_000)

# 5) Get radius
# (returns radius in Å for atom or ion)
crystal.get_radius("Li")                # returns radius for Li atom
crystal.get_radius("Ca", 2, "VIII") # returns ionic radius for Ca2+ ion (VIII CN)

# 6) Get neighbors
# (returns list of neighbors for a given position and atom type)
# tolerance allows for some flexibility in distance matching
# system can be "int", "frac", or "cart" for integer, fractional, or cartesian coordinates
# ball=True returns all neighbors within cutoff + tolerance
# ball=False returns neighbors within [cutoff - tolerance, cutoff + tolerance]
crystal.get_neighbors(0, 0, 0, cutoff=3.0, tolerance=0.1, system="int", pos_rounding="int", ball=True)

# 7) Get neighbors with specific atom type
# (returns list of neighbors for a given position and atom type)
# this gets the cnf variable IDs for neighbors of Li at (0,0,0) within a cutoff of 3.0 Å
crystal.get_neighbors(0, 0, 0, cutoff=3.0, tolerance=0.1, system="int", pos_rounding="int", ball=True, atom_id=crystal.atom_id("Li"))

#--- dictionaries ----#
# variable dictionary
crystal.var_dict[0] # returns the atom and position for variable ID 0
# reverse variable dictionary
crystal.reverse_var_dict[1,1,1 , "Li"] # returns the variable ID for Li at position (1,1,1)

# ion dictionary
crystal.ion_dict[6] # returns the variable ID for Ca2+ ion (VIII CN)
# reverse ion dictionary
crystal.reverse_ion_dict["Ca", 2, "VIII"] # returns the variable ID for Ca2+ ion (VIII CN)


#---- How to use constraints ----#

# 1) Fix atoms (optional)
crystal.force_atom_at_position(0, 0, 0, crystal.atom_id("Li"))
crystal.force_atom_at_position(1, 1, 1, crystal.atom_id("C"))

# 2) Count bounds
crystal.bound_atom(0, 10, 20)   # Li
crystal.bound_atom(2, 10, 20)   # C
crystal.bound_atom(4, 7, 20)    # Ca2+ (VIII)

# 3) Fill all sites (optional)
crystal.fill_unit_cell()

# 4) Geometry constraints (optional)
li_id = crystal.atom_id("Li")
ca_id = crystal.atom_id("Ca", 2, "VIII")
ru_id = crystal.atom_id("Ru", 8, "IV")

crystal.isolate(ru_id, cutoff=2.5, tolerance=0.1)
crystal.isolate_from_types(li_id, [ca_id], cutoff=2.0, tolerance=0.1)
crystal.enforce_closest_dist(li_id, ca_id, min_dist=0.5)

# 5) Solve
solutions = crystal.solve_multiple(n_solutions=10)

# 6) Decode + export
if solutions:
    print("\n✅ Solutions found.")
    for i, sol in enumerate(solutions, start=1):
        decoded = crystal.decode_solution(sol, system_output="int")
        print(f"\nSolution {i}:")
        for entry in decoded:
            print("  ", entry)
        crystal.export_to_CIF(sol, f"test_{i}.cif")
        print(f"   → wrote test_{i}.cif")
else:
    print("\n❌ No solution found.")

# Quick helper reference:
#   atom_id(symbol[, charge, cn])
#   inverse_id(atom_id)
#   get_radius(...)
#   force_atom_at_position(...)
#   forbid_atom_at_position(...)
#   bound_atom(...)
#   get_neighbors(...)