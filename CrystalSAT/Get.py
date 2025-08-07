
import json
from mendeleev import element

class GetMixin:

    def get_positions(self, atom_id):
        """
        Grabs all variable IDs associated with a given atom type.
        :param atom_id: ID of the atom type
        :return: array of variable IDs
        """
        return [
            self.encode_var(x, y, z, atom_id)
            for x,y,z in self.positions
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
        return [self.encode_var(x_int, y_int, z_int, t) for t in range(self.lower,self.k)]

    @staticmethod
    def get_ions(symbol):
        """
        Gets all possible ions for a given element symbol based on Shannon radii database.
        :param symbol:
        :return: ion array of tuples (charge, coordination number)
        """

        with open("shannon-radii.json") as fp:
            shannon_data = json.load(fp)

        if symbol not in shannon_data:
            raise KeyError(f"Symbol {symbol} not found in shannon radii database")

        ions = []

        for charge_str, cn_dict in shannon_data[symbol].items():
            try:

                charge = int(charge_str)

            except ValueError:
                continue

            for cn_str in cn_dict.keys():

                ions.append(( charge, cn_str))

        return ions


    def get_radius(self, symbol, charge= None, cn = None):
        """
        Gets radius of atom or ion based on its symbol, charge, and coordination number.
        :param symbol: str element symbol (e.g., "O", "Na", "Fe")
        :param charge: int charge of the ion (for ions)
        :param cn: coordination number (for ions)
        :return: float radius in Angstroms
        """

        if charge is None and cn is None:
            el = element(symbol)
            r = el.vdw_radius

            if r is None:
                r = el.covalent_radius
                if r is None:
                    raise ValueError(f"No radius found for element {symbol}.")

            r_angstrom = r / 100
            return r_angstrom

        elif charge is not None and cn is not None:

            key = ( symbol, charge, cn)
            if key not in self.reverse_ion_dict:
                raise ValueError(f"Symbol {symbol} is not a valid element or does not have a defined radius.Error : {key}")

            with open("CrystalSATv2/shannon-radii.json") as fp:
                shannon_data = json.load(fp)

            try:
                r_dict = shannon_data[symbol][str(charge)][cn]
                r_ionic = r_dict.get("r_ionic")
                if r_ionic is None:
                    raise ValueError(f"No Shannon ionic radius found for {key}. ")
                return r_ionic
            except KeyError as e:
                raise ValueError(f"No Shannon ionic radius found for {key}. Error: {e}")

        else:
            raise ValueError("Either both charge and cn must be specified (ions), or neither (atoms) .")


    def get_max_radius(self):
        """
        Gets the maximum radius of all particles in the unit cell.
        :return: float maximum radius in Angstroms
        """

        max_radius = 0.0
        for i in range(self.lower,self.k):
            particle = self.inverse_id(i)
            if isinstance(particle, tuple):
                r = self.get_radius(*particle)
            else:
                r = self.get_radius(particle)
            if r > max_radius:
                max_radius = r
        return max_radius






