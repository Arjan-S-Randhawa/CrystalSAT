from pysat.card import CardEnc

class CardinalityMixin:

    def bound_atom(self, atom_id, min_count = None, max_count = None):

        available  = self.grab_available_positions(atom_id)
        forced_count = len(self.grab_forced(atom_id))
        if min_count is not None:
            adjusted_min_count = min_count
            if adjusted_min_count > forced_count + len(available):
                raise ValueError(f"min_count {min_count} is greater than the number of available positions {len(available)} + forced atoms {forced_count} for atom type {atom_id}.")

            atleast = CardEnc.atleast(lits = available, bound=adjusted_min_count, encoding = 1, vpool= self.vpool)

            self.cnf.extend(atleast.clauses)
            self.vpool.start_from = self.vpool.id

        if max_count is not None:
            adjusted_max_count = max_count - forced_count
            if adjusted_max_count < 0:
                raise ValueError(f"max_count {max_count} is less than the number of forced atoms {forced_count} for atom type {atom_id}.")
            else:
                atmost = CardEnc.atmost(lits = available, bound=adjusted_max_count, encoding = 1, vpool = self.vpool)
                self.cnf.extend(atmost.clauses)
                self.vpool.start_from = self.vpool.id

