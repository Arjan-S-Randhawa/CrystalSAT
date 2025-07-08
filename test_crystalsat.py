
# Add the parent directory (Main) to Python's module search path


from CrystalSAT import CrystalSAT

def test_basic_crystal():

    # Create 2x2x2 grid with Lithium and Oxygen
    allowed_atoms = ["Li","O","He","H"]
    model = CrystalSAT(n=4, allowed=[])

    # Initialize constraints: at most 1 atom per site
    model.initialise()

    # Force Fe at (1, 1, 1)
    # problem at 3
    model.force_atom_at_position(
        x=1, y=1, z=1,
        atom_id=model.atom_id("Fe"),
        pos_rounding="int",
        frac=False
    )

    # Force Fe at (0, 0, 0)
    model.force_atom_at_position(
        x=0, y=0, z=0,
        atom_id=model.atom_id("Fe"),
        pos_rounding="int",
        frac=False
    )
    # Force Fe at (2,2,2)
    model.force_atom_at_position(
        x = 2, y = 2, z =2,
        atom_id= model.atom_id("Fe"),
        pos_rounding="int",
        frac=False
    )


    model.at_most(amount = 2, atom_id= model.atom_id("Fe"))









    # Solve
    solution = model.solve()
    if solution:
        print("✅ SATISFIABLE: Solution found.")
        decoded = model.decode_solution(solution, frac=False)
        for entry in decoded:
            print(entry)  # (x, y, z, atom_symbol, truth_value)

    else:
        print("❌ UNSATISFIABLE: No solution found.")


if __name__ == "__main__":
    test_basic_crystal()
