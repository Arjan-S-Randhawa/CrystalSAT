class GrabMixin:

    def grab_forced(self,atom_id):

        """
        :param atom_id:
        :return: list of forced true variables for a specific atom type

        returns all the vars that have been forced true for a specific type
        helper function for global cardinality constraints
        """
        forced_true_vars = [
            clause[0] for clause in self.cnf.clauses
            if len(clause) == 1 and clause[0] > 0 and (clause[0] - 1) % self.k == atom_id

        ]
        return forced_true_vars

    # Count Mixin
    def grab_available_positions(self, atom_id):
        """
        Gets all available positions for a specific atom type.
        Does this by removing occupied positions across all atom types from all possible positions.
        :param atom_id:
        :return: list of available positions for a specific atom type
        """

        forced_positions = set()
        for t in range(self.k):
            for var in self.grab_forced(t):

                x,y,z,_ = self.decode_var(var)

                forced_positions.add((x,y,z))

        occupied = set(
            self.encode_var(x,y,z,atom_id)
            for (x,y,z) in forced_positions
        )

        all_vars = set(self.get_positions(atom_id))

        available_positions = all_vars - occupied

        return list(available_positions)

