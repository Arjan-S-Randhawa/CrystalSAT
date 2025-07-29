
class ConstraintsMixin:

    def initialise(self):
        """
        Initialises the CNF with basic constraints.
        Ensures no two atoms of the same type can occupy the same position.
        :return:
        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    types_for_this_pos = self.get_types(x,y,z, system="int", pos_rounding="int")
                    for i in range(self.k):
                        for j in range(i + 1, self.k):
                            self.cnf.append([-types_for_this_pos[i], -types_for_this_pos[j]])



    def fill_unit_cell(self):

        """
        Forces CNF to fill all positions in the unit cell with atom types.
        :return:

        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    types = self.get_types(x, y, z)
                    self.cnf.append(types)



    def force_atom_at_position(self, x,y,z,atom_id,system ="int", pos_rounding= "int"):
        """
        :param x:
        :param y:
        :param z:
        :param atom_id:
        :param system:
        :param pos_rounding:
        :return:
        """
        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        self.cnf.append([self.encode_var(x_int, y_int, z_int, atom_id)])


    def forbid_atom_at_position(self, x,y,z,atom_id,system = "int",pos_rounding = "int"):
        """
        :param x:
        :param y:
        :param z:
        :param atom_id:
        :param system:
        :param pos_rounding:
        :return:
        """
        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        self.cnf.append([-self.encode_var(x_int,y_int,z_int,atom_id)])


    def require_one_of_types_at_position(self, x,y,z,atom_ids, system = "int", pos_rounding = "int"):
        """
        :param x:
        :param y:
        :param z:
        :param atom_ids: array []
        :param system:
        :param pos_rounding:
        :return:
        """

        x_int, y_int, z_int = self.to_int(x,y,z,system=system,pos_rounding=pos_rounding)
        vars = [self.encode_var(x_int,y_int,z_int,atom_id) for atom_id in atom_ids]
        self.cnf.append(vars)


    def forbid_types_at_position(self, x, y, z, atom_ids, system = "int", pos_rounding = "int"):
        """
        Forbids all specified atom types at a given position (x, y, z).
        :param x:
        :param y:
        :param z:
        :param atom_ids:
        :param system:
        :param pos_rounding:
        :return:
        """
        for atom_id in atom_ids:
            self.forbid_atom_at_position(x,y,z,atom_id, system = system, pos_rounding = pos_rounding)
