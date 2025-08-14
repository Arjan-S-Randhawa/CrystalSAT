from CrystalSAT import *
# Example: PbTiO3 in Pm-3m without packing

# 1) Build the cell (cubic, Pm-3m)
a = b = c = 3.97
alpha = beta = gamma = 90.0

# specify the allowed ion types, their charges, and coordination numbers
allowed = [
    ('Pb',  2, 'XII'),  # A site
    ('Ti',  4, 'VI'),   # B site
    ('O',  -2, 'II')    # oxygen
]
#create
crystal = CrystalSAT(n_x=2, n_y=2, n_z=2,
                 a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma,
                 allowed=allowed)

# initialise the model, turn the sat solver on
# pack is set to false for this
crystal.initialise(pack=False)

# 3) Get type IDs
# sat solvers need atom ids to encode constraints, this is how crystalSAT knows what atom you are referring to
pb = crystal.atom_id('Pb', 2, 'XII')
ti = crystal.atom_id('Ti', 4, 'VI')
ox = crystal.atom_id('O', -2, 'II')

# ensures no Pb can be placed next to each other
crystal.isolate_from_itself(pb,a/2,0.1,ball = True)

# this forces titanium to be surrounded by oxygen atoms
# a/2 is the distance specified, 0.1 is the tolerance

crystal.surround_type(ti,ox,a/2,0.1)

# 5) Exact counts
crystal.bound_atom(pb, min_count=1, max_count=1)
crystal.bound_atom(ti, min_count=1, max_count=1)
crystal.bound_atom(ox, min_count=3, max_count=3)


# 7) Solve and export
# enumarates over 10,0000 solutions and removes symmetrical solutioms
# this includes PbTiO3 in Pm-3m
# ends up being only 2 unique unit cells
unique_solutions,uniques, groups= crystal.solve_unique(n_solutions =10000)
if unique_solutions:
    print("Solution found!")
    print("number of uniques", len(uniques))
    print("groups", len(groups))
    for i,sol in enumerate(unique_solutions):
        print("Solution:", i+1)
        decode = crystal.decode_solution(sol, system_output='frac')
        for decoded in decode:
            print(decoded)

        crystal.export_to_CIF(sol, filename=f'PbTiO3_pm3m_variation_{i+1}.cif')

else:
    print("No solution found.")
