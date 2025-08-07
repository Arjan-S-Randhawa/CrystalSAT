from ase.io import write
from ase import Atoms

class SolveAndExportMixin:


    def solve(self, solver_name="glucose3"):
        """
        Solves the SAT problem using the specified solver and returns a model if satisfiable.
        :param solver_name:
        :return: model if satisfiable, None if unsatisfiable
        """
        from pysat.solvers import Solver
        with Solver(name=solver_name, bootstrap_with=self.cnf.clauses) as solver:
            is_sat = solver.solve()
            if is_sat:
                return solver.get_model()
            else:
                return None

    def solve_multiple(self, solver_name="glucose3", n_solutions=1):
        """
        Solves the SAT problem and returns multiple solutions.
        :param solver_name: Name of the SAT solver to use
        :param n_solutions: Number of solutions to find
        :return: List of solutions, where each solution is a list of integers representing the model
         """
        from pysat.solvers import Solver

        solutions = []
        cnf = [clause[:] for clause in self.cnf.clauses]  # Make a copy

        with Solver(name=solver_name, bootstrap_with=cnf) as solver:
            for _ in range(n_solutions):
                is_sat = solver.solve()
                if not is_sat:
                    break
                model = solver.get_model()
                solutions.append(model)

                # Create a blocking clause to prevent this exact solution from repeating
                blocking_clause = [-lit if lit > 0 else -lit for lit in model if abs(lit) in self.var_dict]
                solver.add_clause(blocking_clause)

        return solutions

    def decode_solution(self, solution, system_output="int"):
        """
        Decodes a solution from the SAT solver into a list of variables.
        :param solution: List of integers representing the SAT solution
        :return: list of true variables in the format (x, y, z, atom_symbol, truth_value)
        """
        max_original = self.n_x * self.n_y * self.n_z * self.k
        variables = []

        for encoded_var in solution:

            if 0 < encoded_var <= max_original and encoded_var in self.var_dict:
                x, y, z, atom_symbol = self.var_dict[encoded_var]

                if system_output == "frac":
                    x_s, y_s, z_s = self.to_frac(x, y, z, system="int", pos_rounding="int")
                elif system_output == "cart":
                    x_s, y_s, z_s = self.to_cart(x, y, z, system="int", pos_rounding="int")
                elif system_output == "int":
                    x_s, y_s, z_s = x, y, z

                variables.append((x_s, y_s, z_s, atom_symbol))

        return variables


    def export_to_ase(self, solution):
        """
        Exports the solution to an ASE Atoms object.
        :param solution: SAT solution to export
        :return: ASE Atoms object
        """

        decoded = self.decode_solution(solution, system_output="cart")
        symbols = []
        charges = []
        cart_positions = []

        for x, y, z, atom_symbol in decoded:

            if isinstance(atom_symbol, tuple):
                symbol,charge, cn = atom_symbol
                ase_symbol = symbol
                ase_charge = charge

            else:
                ase_symbol = atom_symbol
                ase_charge = 0

            symbols.append(ase_symbol)
            charges.append(ase_charge)
            cart_positions.append([x, y, z])

        atoms = Atoms( symbols = symbols, positions = cart_positions , cell=self.cell, pbc=True)
        return atoms


    def export_to_CIF(self, solution, filename="output.cif"):
        """
        Exports a solved SAT model to a CIF file.
        :param solution: SAT solution to export
        :param filename: Name of the output CIF file
        :return: None
        """
        atoms = self.export_to_ase(solution)
        write(filename, atoms, format='cif')
        print(f"CIF file is saved to {filename}")


