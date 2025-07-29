
import periodictable

class EncodingMixin:

    @staticmethod
    def atomic_num(symbol):
        """
        :param symbol: str
        :return: int
        converts symbol to atomic number
        """
        return periodictable.elements.symbol(symbol).number

    # encode mixin
    def atom_id(self, symbol):
        """
        :param symbol: str
        :return: int

        Assigns each atom type a unique int.
        """
        if self.use_allowed:
            if symbol not in self.allowed:
                raise ValueError(f"Symbol {symbol} is not in the allowed list: {self.allowed}")
            else:
                return self.allowed.index(symbol)
        else:
            return self.atomic_num(symbol)

    # encode mixin
    def inverse_id(self, atom_id):
        """
        :param atom_id: integer ID of the atom type
        :return: string symbol of the atom type
        """
        if self.use_allowed:
            if not (0 <= atom_id < len(self.allowed)):

                raise ValueError(f"Atom ID {atom_id} is out of bounds for allowed atom types.")
            else:
                return self.allowed[atom_id]
        else:
            if not (1 <= atom_id <= 118):
                raise ValueError(f"Atom ID {atom_id} is out of bounds for the periodic table.")
            else:
                return periodictable.elements[atom_id].symbol

    # encode mixin
    def encode_var(self, x, y, z, atom_id):
        """
        Encodes a 3D position (x, y, z) and an atom type (k) into a unique integer.
        :param x: x-coordinate (integer)
        :param y: y-coordinate (integer)
        :param z: z-coordinate (integer)
        :param atom_id: atom type (integer ID)
        :return:Integer ID representing the variable
        """
        var_id = (
                         x * (self.n_y * self.n_z * self.k) +
                         y * (self.n_z * self.k) +
                         z * self.k +
                         atom_id
                 ) + 1
        return int(var_id)

    # encode mixin
    def decode_var(self, var_id):
        idx = var_id - 1
        atom_id = idx % self.k
        idx //= self.k
        z = idx % self.n_z
        idx //= self.n_z
        y = idx % self.n_y
        idx //= self.n_y
        x = idx
        return x, y, z, atom_id


    def populate_var_dict(self):
        """
        Populates the variable dictionary with all possible variables.
        :return: None
        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    for k in range(self.k):
                        var_id = self.encode_var(x, y, z, k)
                        self.var_dict[var_id] = (x, y, z, self.inverse_id(k))

