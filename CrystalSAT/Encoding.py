
import periodictable
import json

class EncodingMixin:

    @staticmethod
    def atomic_num(symbol):
        """
        Converts an element symbol to its atomic number.
        :param symbol: Element symbol (e.g., 'H', 'O', 'Fe')
        :return: Atomic number of the element (int)
        """
        return periodictable.elements.symbol(symbol).number


    def populate_ion_dict(self):
        """
        Populates the ion dictionary with elements and their charges from the Shannon radii data.
        :return:
        """
        idx = 119
        with open("CrystalSATv2/shannon-radii.json") as fp:
            shannon_data = json.load(fp)

        for symbol, charges in shannon_data.items():
            for charge_str, cn_dict in charges.items():
                try:
                    charge = int(charge_str)
                except ValueError:
                    continue
                for cn_str in cn_dict.keys():
                    self.ion_dict[idx] = (symbol, charge, cn_str)
                    self.reverse_ion_dict[(symbol, charge, cn_str)] = idx

                    idx+=1


    def atom_id(self,symbol,charge = 0,cn = None):
        """
        Converts an atom symbol, charge, and coordination number to a unique integer ID.
        :param symbol:
        :param charge:
        :param cn:
        :return: int ID of the atom type
        """

        if self.use_allowed:
            # atom logic
            if charge == 0 and cn is None:
                if symbol in self.allowed:
                    return self.allowed.index(symbol)
                else:
                    raise ValueError(f"Atom symbol {symbol} is not in the allowed list.")

            elif charge != 0 and cn is not None:

                for i, item in enumerate(self.allowed):
                    if isinstance(item,tuple) and item==(symbol,charge,cn):
                        return i

                raise ValueError(f"Atom symbol {symbol} with charge {charge} and coordination number {cn} is not in the allowed list.")

            else:
                raise ValueError("Ions require both charge and coordination number to be specified.")

        else:
            if charge == 0 and cn is None:
                return self.atomic_num(symbol)

            elif charge != 0 and cn is not None:

                key = (symbol,charge,cn)
                for idx, val in self.ion_dict.items():
                    if val == key:
                        return idx
                raise ValueError(f"Atom symbol {symbol} with charge {charge} and coordination number {cn} is not in the ion dictionary.")
            else:
                raise ValueError("Ions require both charge and coordination number to be specified.")



    def inverse_id(self, atom_id):
        """
        Converts an atom type ID back to its symbol or ion representation.
        :param atom_id: integer ID of the atom type
        :return: string symbol of the atom type or tuple for ion
        """
        if self.use_allowed:
            if not (0 <= atom_id < len(self.allowed)):
                raise ValueError(f"Atom ID {atom_id} is out of bounds for allowed atom types.")
            else:
                return self.allowed[atom_id]
        else:
            if not (1 <= atom_id <= self.k):
                raise ValueError(f"Atom ID {atom_id} is out of bounds for the periodic table.")
            else:
                if (1<= atom_id <= 118):

                    return periodictable.elements[atom_id].symbol

                elif (119<= atom_id <= self.k):

                    return self.ion_dict[atom_id]

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

    def populate_var_dict(self):
        """
        Populates the variable dictionary with all possible variables.
        :return: None
        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    for k in range(self.lower,self.k):
                        var_id = self.encode_var(x, y, z, k)
                        self.var_dict[var_id] = (x, y, z, self.inverse_id(k))
                        self.reverse_var_dict[(x, y, z, self.inverse_id(k))] = var_id



