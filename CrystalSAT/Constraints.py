
from itertools import combinations, combinations_with_replacement

class ConstraintsMixin:


    def initialise(self, pack = True):
        """
        Initialises the CNF with basic constraints.
        Ensures no two atoms can occupy the same position.
        Ensures that atoms do not overlap based on their radii (if pack is True).
        :return:
        """

        for x,y,z in self.positions:
            # needs fixing
            types_for_this_pos = self.get_types(x, y, z, system="int", pos_rounding="int")
            # No two atoms of the different type can occupy the same position
            for i,j in combinations(range(self.lower,self.k), 2):

                self.cnf.append([-types_for_this_pos[i-self.lower], -types_for_this_pos[j-self.lower]])

        # Sphere packing constraints: Forbids overlapping atoms based on their radii.
        # Uses ionic radii for ions and vdw/covalent radii for atoms

        if pack:
            max_radius = self.get_max_radius()
            radii = [self.get_radius(*self.inverse_id(i)) if isinstance(self.inverse_id(i), tuple) else self.get_radius(
                self.inverse_id(i)) for i in range(self.lower,self.k)]
            cutoff = max_radius * 2.0

            for (x1, y1, z1), (x2, y2, z2) in combinations(self.positions, 2):
                idx1 = x1 * (self.n_y * self.n_z) + y1 * self.n_z + z1
                idx2 = x2 * (self.n_y * self.n_z) + y2 * self.n_z + z2

                dist = self.grid.get_distance(idx1, idx2, mic=True)
                if dist > cutoff:
                    continue

                for i, j in combinations_with_replacement(range(self.lower,self.k), 2):
                    rad_1 = radii[i - self.lower]
                    rad_2 = radii[j - self.lower]

                    if dist < rad_1 + rad_2:
                        self.cnf.append([-self.encode_var(x1, y1, z1, i),
                                         -self.encode_var(x2, y2, z2, j)])
                        if i != j:
                            self.cnf.append([-self.encode_var(x1, y1, z1, j),
                                            -self.encode_var(x2, y2, z2, i)])


    def fill_unit_cell(self):

        """
        Forces solver to fill all positions in the unit cell with atom types.
        :return:
        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    types = self.get_types(x, y, z)
                    self.cnf.append(types)



    def force_atom_at_position(self, x,y,z,atom_id,system ="int", pos_rounding= "int"):
        """
        Forces a specific atom type to occupy a given position (x, y, z).

        :param x: x-coordinate of the position
        :param y: y-coordinate of the position
        :param z: z-coordinate of the position
        :param atom_id: ID of the atom type to force
        :param system: specifies the coordinate system to use
        :param pos_rounding: when typing in invalid cartesian/ fractional coordinates, this parameter specifies how to round them
        :return:
        """
        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        self.cnf.append([self.encode_var(x_int, y_int, z_int, atom_id)])


    def forbid_atom_at_position(self, x,y,z,atom_id,system = "int",pos_rounding = "int"):
        """
        Forbids a specific atom type from occupying a given position (x, y, z).

        :param x: x-coordinate of the position
        :param y: y-coordinate of the position
        :param z: z-coordinate of the position
        :param atom_id: ID of the atom type to forbid
        :param system: specifies the coordinate system to use
        :param pos_rounding: when typing in invalid cartesian/ fractional coordinates, this parameter specifies how to round them
        :return:
        """
        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        self.cnf.append([-self.encode_var(x_int,y_int,z_int,atom_id)])


    def require_one_of_types_at_position(self, x,y,z,atom_ids, system = "int", pos_rounding = "int"):
        """
        Requires at least one of the specified atom types to occupy a given position (x, y, z).

        :param x: x-coordinate of the position
        :param y: y-coordinate of the position
        :param z: z-coordinate of the position
        :param atom_ids: array []
        :param system: specifies the coordinate system to use
        :param pos_rounding: when typing in invalid cartesian/ fractional coordinates, this parameter specifies how to round them
        :return:
        """

        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        vars = [self.encode_var(x_int,y_int,z_int,atom_id) for atom_id in atom_ids]
        self.cnf.append(vars)


    def forbid_types_at_position(self, x, y, z, atom_ids, system = "int", pos_rounding = "int"):
        """
        Forbids all specified atom types at a given position (x, y, z).

        :param x: x-coordinate of the position
        :param y: y-coordinate of the position
        :param z: z-coordinate of the position
        :param atom_ids: list of atom IDs to forbid at the position
        :param system: specifies the coordinate system to use
        :param pos_rounding: when typing in invalid cartesian/ fractional coordinates, this parameter specifies how to round them
        :return:
        """
        for atom_id in atom_ids:
            self.forbid_atom_at_position(x,y,z,atom_id, system = system, pos_rounding = pos_rounding)
