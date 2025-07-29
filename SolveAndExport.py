from ase.io import write
from ase import Atoms

class SolveAndExportMixin:


    def solve(self, solver_name="glucose3"):
        from pysat.solvers import Solver

        with Solver(name=solver_name, bootstrap_with=self.cnf.clauses) as solver:
            is_sat = solver.solve()
            if is_sat:
                return solver.get_model()
            else:
                return None


    def decode_solution(self, solution, system_output="int"):
        """
        Decodes a solution from the SAT solver into a list of variables.
        :param solution:
        :return: list of true variables in the format (x, y, z, atom_symbol, truth_value)
        """
        max_original = self.n_x * self.n_y * self.n_z * self.k
        variables = []

        for encoded_var in solution:

            if 0 < encoded_var <= max_original:
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
        :param self:
        :param solution:
        :return: ASE Atoms object
         """

        decoded = self.decode_solution(solution, system_output="cart")
        symbols = []
        cart_positions = []

        for x, y, z, atom_symbol in decoded:
            symbols.append(atom_symbol)
            cart_positions.append([x, y, z])

        atoms = Atoms(symbols=symbols, positions=cart_positions, cell=self.cell, pbc=True)
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


