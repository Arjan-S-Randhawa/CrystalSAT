from math import *
from numpy import *


class CoordinateMixin:

    def to_int(self, x, y, z, system = "int", pos_rounding="int"):
        """
        Converts fractional or cartesian coordinates to integer coordinates.
        :param x: fractional/ cartesian x-coordinate
        :param y: fractional/cartesian y-coordinate
        :param z: fractional z-coordinate
        :param system: coordinate system to use ("frac" or "cart")
        :param pos_rounding: how to round fractional positions ("floor", "ceil", or "round")
        :return: int tuple (x, y, z)

        converts fractional coordinates to cartesian coordinates.
        """
        x_int = y_int = z_int = None

        if system == "frac":

            max_frac_x = (self.n_x - 1) / self.n_x
            max_frac_y = (self.n_y - 1) / self.n_y
            max_frac_z = (self.n_z - 1) / self.n_z

            if not 0 <=x < 1 or not 0<= y < 1 or not 0<= z < 1:

                raise ValueError(f"Fractional coordinates ({x}, {y}, {z}) must be between 0 and 1")

            else:

                x = min(x, max_frac_x)
                y = min(y, max_frac_y)
                z = min(z, max_frac_z)

                if pos_rounding == "floor":

                    x_int = int(x * self.n_x)
                    y_int = int(y * self.n_y)
                    z_int = int(z * self.n_z)

                elif pos_rounding == "ceil":

                    x_int = int(ceil(x * self.n_x))
                    y_int = int(ceil(y * self.n_y))
                    z_int = int(ceil(z * self.n_z))

                elif pos_rounding == "round":

                    x_int = int(round(x * self.n_x))
                    y_int = int(round(y * self.n_y))
                    z_int = int(round(z * self.n_z))
                else:

                    raise ValueError("Unsupported rounding type")

        elif system == "cart":

            if not (0 <= x < self.a and 0 <= y < self.b and 0 <= z < self.c):
                raise ValueError(
                    f"Cartesian coordinates ({x}, {y}, {z}) must be within the unit cell dimensions ({self.a}, {self.b}, {self.c})")

            else:

                dx = self.a / self.n_x
                dy = self.b / self.n_y
                dz = self.c / self.n_z
                max_cart_x = (self.n_x - 1) * dx
                max_cart_y = (self.n_y - 1) * dy
                max_cart_z = (self.n_z - 1) * dz

                x = min(x, max_cart_x)
                y = min(y, max_cart_y)
                z = min(z, max_cart_z)

                if pos_rounding == "floor":
                    x_int = int(x/ dx)
                    y_int = int(y / dy)
                    z_int = int(z / dz)

                elif pos_rounding == "ceil":
                    x_int = int(ceil(x / dx))
                    y_int = int(ceil(y / dy))
                    z_int = int(ceil(z / dz))

                elif pos_rounding == "round":
                    x_int = int(round(x / dx))
                    y_int = int(round(y / dy))
                    z_int = int(round(z / dz))

        elif system == "int" or pos_rounding == "int":

            if not isinstance(x, int) or not isinstance(y, int) or not isinstance(z, int):
                raise ValueError("Integer coordinates must be provided as integers")
            else:
                if not (0 <= x < self.n_x and 0 <= y < self.n_y and 0 <= z < self.n_z):
                    raise ValueError(f"Integer coordinates ({x}, {y}, {z}) are not within the grid dimensions ({self.n_x}, {self.n_y}, {self.n_z})")
                else:
                    x_int, y_int, z_int = x, y, z

        return x_int, y_int, z_int

    def to_frac(self, x, y, z, system = "int" , pos_rounding = "int"):
        """
        Converts integer or cartesian coordinates to fractional coordinates.

        :param x: x-coordinate (integer or cartesian)
        :param y: y-coordinate (integer or cartesian)
        :param z: z-coordinate (integer or cartesian)
        :param system: coordinate system to use ("int", "frac", or "cart")
        :param pos_rounding: how to round fractional positions ("floor", "ceil", or "round")
        :return: tuple coordinates
        """

        if system == "frac":

            x_int, y_int, z_int = self.to_int(x, y, z, system=system, pos_rounding=pos_rounding)
            x_frac = x_int / self.n_x
            y_frac = y_int / self.n_y
            z_frac = z_int / self.n_z

        else:

            x_int, y_int, z_int = self.to_int(x, y, z, system, pos_rounding)
            x_frac = x_int / self.n_x
            y_frac = y_int / self.n_y
            z_frac = z_int / self.n_z


        return x_frac, y_frac, z_frac

    def to_cart(self,x,y,z, system="int",pos_rounding= "int"):

        """
        Converts fractional or integer coordinates to cartesian coordinates.

        :param x: x-coordinate (fractional or integer)
        :param y: y-coordinate (fractional or integer)
        :param z: z-coordinate (fractional or integer)
        :param system: coordinate system to use ("int", "frac", or "cart")
        :param pos_rounding: how to round fractional positions ("floor", "ceil", or "round")
        :return: x,y,z in cartesian coordinates
        """

        dx = self.a / self.n_x
        dy = self.b / self.n_y
        dz = self.c / self.n_z

        x_int, y_int, z_int = self.to_int(x, y, z, system = system , pos_rounding = pos_rounding)
        x_cart = x_int * dx
        y_cart = y_int * dy
        z_cart = z_int * dz

        return x_cart, y_cart, z_cart
