
class GetMixin:

    def get_positions(self, atom_id):
        """
        Grabs all variable IDs associated with a given atom type.
        :param atom_id:
        :return: array of variable IDs
        """
        return [
            self.encode_var(x, y, z, atom_id)
            for x in range(self.n_x)
            for y in range(self.n_y)
            for z in range(self.n_z)
        ]

    def get_types(self, x, y, z, system = "int", pos_rounding = "int"):
        """
        Grabs all variable IDs associated with a given position (x, y, z).
        :param x: x-coordinate (integer)
        :param y: y-coordinate (integer)
        :param z: z-coordinate (integer)
        :param system: coordinate system to use ("int", "frac", or "cart")
        :param pos_rounding: how to round fractional positions ("floor", "ceil", or "round")
        :return: array of variable IDs
        """
        x_int, y_int, z_int = self.to_int(x, y, z, system, pos_rounding)


        return [self.encode_var(x_int, y_int, z_int, t) for t in range(self.k)]




