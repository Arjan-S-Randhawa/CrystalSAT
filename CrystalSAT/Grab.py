class GrabMixin:

    def grab_forced(self,atom_id):

        """
        Gets all forced true variables for a specific atom type.
        :param atom_id:
        :return: list of forced true variables for a specific atom type
        """
        forced_true_vars = [
            clause[0] for clause in self.cnf.clauses
            if len(clause) == 1 and clause[0] > 0 and (clause[0] - 1) % self.k == atom_id

        ]
        return forced_true_vars

    # Count Mixin
    def grab_available_positions(self, atom_id):
        """
        Gets all available positions for a specific atom type, excluding those that are already occupied by forced atoms.
        :param atom_id: ID of the atom type to check for available positions
        :return: list of available positions for a specific atom type
        """
        forced_positions = set()
        for t in range(self.lower,self.k):
            for var in self.grab_forced(t):

                x,y,z,_ = self.var_dict[var]

                forced_positions.add((x,y,z))

        occupied = set(
            self.encode_var(x,y,z,atom_id)
            for (x,y,z) in forced_positions
        )

        all_vars = set(self.get_positions(atom_id))

        available_positions = all_vars - occupied

        return list(available_positions)

