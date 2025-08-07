
from ase.neighborlist import NeighborList


class NeighborAndDistancesMixin:

    def get_distance(self, x1, y1, z1, x2, y2, z2):

        """
        Calculates the distance between two points in the grid defined
        by their coordinates (x1, y1, z1) and (x2, y2, z2).

        :param x1: x-coordinate of the first point
        :param y1: y-coordinate of the first point
        :param z1: z-coordinate of the first point
        :param x2: x-coordinate of the second point
        :param y2: y-coordinate of the second point
        :param z2: z-coordinate of the second point
        :return: distance between the two points in Angstroms
        """

        idx1 = x1 * (self.n_y * self.n_z) + y1 * self.n_z + z1
        idx2 = x2 * (self.n_y * self.n_z) + y2 * self.n_z + z2
        return self.grid.get_distance(idx1, idx2, mic=True)


    def get_neighbors(self, x, y, z, cutoff , tolerance, system , pos_rounding, debug = True , ball = True):
        """
        Finds all neighbors of a given position (x, y, z) within a specified cutoff distance (Å).
        Input can be integer, fractional, or cartesian; output is a list of (x', y', z') in grid indices.

        :param x: x-coordinate (integer, fractional or cartesian)
        :param y: y-coordinate (integer, fractional or cartesian)
        :param z: z-coordinate (integer, fractional or cartesian)
        :param cutoff: cutoff distance (Å) for neighbors
        :param tolerance: tolerance for distance matching (Å)
        :param system: coordinate system to use ("int", "frac", or "cart")
        :param pos_rounding: how to round fractional/cartesian positions ("floor", "ceil", or "round")
        :param debug: if True, prints debug information
        :param ball: if True, returns all neighbors within cutoff + tolerance; if False, returns neighbors within [cutoff - tolerance, cutoff + tolerance]
        """

        # get corresponding integer coordinates
        x_int, y_int, z_int = self.to_int(x, y, z, system = system , pos_rounding=pos_rounding)

        # flatten to int
        idx = x_int * (self.n_y * self.n_z) + y_int * self.n_z + z_int

        # Builds neighbor list (radius per site)
        radii = [cutoff / 2.0] * len(self.grid)
        nl = NeighborList(radii, self_interaction=False, bothways=True)
        nl.update(self.grid)

        # Query neighbors
        indices, offsets = nl.get_neighbors(idx)

        # Unflatten to (x', y', z')
        neighbors = []
        for i in indices:

            x2 = i // (self.n_y * self.n_z)
            y2 = (i % (self.n_y * self.n_z)) // self.n_z
            z2 = i % self.n_z

            dist = self.grid.get_distance(idx, i, mic=True)

            if ball:
                if 0 <= dist <= cutoff + tolerance:
                    neighbors.append((x2, y2, z2))
                    if debug:
                        print(f"Neighbor: ({x2}, {y2}, {z2}), Distance: {dist:.4f} Å")

            else:
                if  ( cutoff - tolerance)<= dist <= cutoff + tolerance:
                    neighbors.append((x2, y2, z2))
                    if debug:
                        print(f"Neighbor: ({x2}, {y2}, {z2}), Distance: {dist:.4f} Å")

        return neighbors

    # NeighborAndDistances Mixin
    def encode_neighbors(self, x, y, z, atom_id, cutoff, tolerance, system ,pos_rounding ,ball=True):

        """
        Encodes the neighbors of a given position (x, y, z) and atom type (atom_id) within a specified cutoff distance.
        Allows us to return a list of encoded variables for neighbors of x,y,z where they have type atom_id.
        Useful for encoding constraints in SAT problems.

        :param x: x-coordinate (integer, fractional or cartesian)
        :param y: y-coordinate (integer, fractional or cartesian)
        :param z: z-coordinate (integer, fractional or cartesian)
        :param atom_id: Atom type (integer ID)
        :param cutoff: Cutoff distance (Å) for neighbors
        :param tolerance: Tolerance for distance matching (Å)
        :param system: Coordinate system to use ("int", "frac", or "cart")
        :param pos_rounding: How to round fractional/cartesian positions ("floor", "ceil", or "round")
        :param ball: If True, returns all neighbors within cutoff + tolerance
        :return: list of encoded neighbor variables for SAT solver
        """
        encoded_neighbors = []
        neighbors = self.get_neighbors(x, y, z, cutoff, tolerance, system=system, pos_rounding=pos_rounding, debug=False, ball=ball)

        for nx, ny, nz in neighbors:
            encoded_neighbors.append(int(self.encode_var(nx, ny, nz, atom_id)))

        return encoded_neighbors



