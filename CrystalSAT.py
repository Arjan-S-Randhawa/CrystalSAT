
"""
CrystalSAT - A Python library for constraint satisfaction modelling in crystal structures
Author: Arjan Randhawa
University: University of Liverpool
Degree: BSc Computer Science and Mathematics
Email: sgarandh@liverpool.ac.uk
This library encodes lattice constraints for use with SAT solvers and integrates with pymatgen for CIF export.
"""

#imports for pysat
from idlelib.autocomplete import FORCE
from tkinter.tix import Select
import periodictable
from pygments.console import *
from pysat.formula import CNF
from itertools import product
from math import *
from pysat.card import *
from numpy import *
from pysat.solvers import Solver
from pymatgen.core import *
import builtins

#TODO
"""
- test current model ( WORKING WELL!!!)
- Expand for more compatible with CIF and 
- possibly implement floor and ceil for spacegroups 
- implement safety features
- fix global constraints 
- change all t/type to atom_type
- switch to maxSAT 



"""

class CrystalSAT:

    def __init__(self, n, allowed):
        """
           creates a CrystalSAT object
           :param n: int
           the number of equally spaced positions allowed 1x1x1 grid
           :param allowed: array[]
           The allowed atom types in the unit cell,
           if left empty all atom types are allowed (118 different atoms)

           """
        self.n = n
        self.allowed = allowed
        self.cnf = CNF()
        self.use_allowed = bool(allowed)

        if self.use_allowed:
            self.k = len(self.allowed)
        else:
            self.k = 118

    # ATOMIC FUNCTIONS
    @staticmethod
    def atomic_num(symbol):
        """
        :param symbol:
        maps pretty symbol to atomic number
        """
        return periodictable.elements.symbol(symbol).number

    def atom_id(self, symbol):
        """
        :param symbol:
        :return: int
        maps pretty symbol to a unique atomic id
        is dependant on if you are using allowed or not
        """
        if self.use_allowed:
            return self.allowed.index(symbol)
        else:
            return self.atomic_num(symbol)- 1  # makes sure it starts from zero

    # CONVERTER FUNCTIONS
    def get_frac(self, x, y, z):
        """
        :param x:
        :param y:
        :param z:
        :return: tuple coordinates

        allows us to turn an integer into fractional coordinates
        no need for rounding
        """
        x_frac = x / self.n
        y_frac = y / self.n
        z_frac = z / self.n
        return x_frac, y_frac, z_frac


    def get_unit(self, x_frac, y_frac, z_frac, pos_rounding = "floor"):
        """
        :param x_frac:
        :param y_frac:
        :param z_frac:
        :param pos_rounding:
        :return: int tuple

        allows us to convert from a fractional to int coordinates
        we can convert any fractional to int
        if fractional coordinates are "not valid positions* we can use round
        to round in accordance to what we want instead
        """
        x = y = z = None
        if pos_rounding == "floor":

            x = int(x_frac * self.n)
            y = int(y_frac * self.n)
            z = int(z_frac * self.n)

        elif pos_rounding == "ceil":

            x = int(ceil(x_frac * self.n))
            y = int(ceil(y_frac * self.n))
            z = int(ceil(z_frac * self.n))

        elif pos_rounding == "round":

            x = int(round(x_frac*self.n))
            y = int(round(y_frac*self.n))
            z = int(round(z_frac*self.n))
        else:

            raise ValueError("Unsupported rounding type")

        return x, y, z

    # ENCODING FUNCTIONS
    def encode_pos(self, x, y, z, pos_rounding , frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param pos_rounding:
        :param frac:
        :return: int

        a bijection of cartesian coordinates to an integer
        we can use fractional coordinates, convert to int system
        if x,y,z is not valid we can round to nearest valid position

        """
        ix = iy = iz = None

        if frac:
            ix,iy,iz = self.get_unit(x,y,z,pos_rounding = pos_rounding)

            if not (0<= x < 1 and 0<= y <1  and 0<= z < 1):
                raise ValueError(f"Fractional coordinates ({x}, {y}, {z}) must be between 0 and 1")

        elif not frac and pos_rounding == "int": # we must pick to use ints
            ix = x
            iy = y
            iz = z
            if not (0<= ix < self.n and 0<= iy < self.n and 0<= iz < self.n):
                raise ValueError(f"Position ({ix},{iy},{iz}) out of bound for grid size n = {self.n} ")
        else:
            raise ValueError("when frac = False, pos_ rounding must be 'int'")
        return ix * self.n * self.n + iy * self.n + iz + 1

    def encode_var(self, x,y,z,atom_id,pos_rounding, frac= True):

        """
        :param x: x-coordinate (fractional or integer, depending on `frac`)
        :param y: y-coordinate
        :param z: z-coordinate
        :param t: atom type (integer ID)
        :param pos_rounding: how to round fractional positions ("floor", "ceil", or "round")
        :param frac: whether coordinates are fractional (True) or integer (False)
        :return: single integer encoding the position and type

        Extends the encode_pos function by including atom type and `k`.
        Maps fractional or integer positions (with type) to a unique integer.

        Example:
        Let n = 4
        Valid fractional positions are:
        0.25, 0.5, 0.75, 1

        We might want to enforce constraints like:
        - A lithium atom must be placed near (0.37, 0.45, 0.76)
        - Or another near (0.67, 0.26, 0.89)

        Since these are *not valid fractional positions*, we generalize:
        They are rounded to the closest valid positions using `pos_rounding`.

        cnf.append([
            -encode_var(0.37, 0.45, 0.76, atom_id("Li"), "round", True),
            -encode_var(0.67, 0.26, 0.89, atom_id("Li"), "round", True)
        ])
        """
        if not (0<= atom_id < self.k):
            raise ValueError(f"Invalid atom tpe {atom_id}. Must be between 0 and {self.k - 1}."
                             f" Please check your allowed atoms list : if you are using atom_id('X') to get a type,"
                             f" ensure 'X' is in your allowed[] atoms.")


        return (self.encode_pos(x,y,z,pos_rounding = pos_rounding,frac = frac) - 1 ) * self.k + atom_id + 1


    # INVERSE FUNCTIONS
    def inverse_id(self, atom_id):
        """
        :param atom_id:
        :return: int

        allows us to get the pretty symbol from atoms id

        """
        if self.use_allowed:
            return self.allowed[atom_id]
        else:
            atomic_number = atom_id +1
            return periodictable.elements[atomic_number].symbol


    def inverse_var(self, var_id ,frac = True):

        """
        :param var_id:
        :param frac:
        :return: 4 tuple
        returns the position and type of a int assigned to it

        """
        var_id -= 1
        pos_id = var_id // self.k
        atom_id = var_id % self.k
        x = pos_id // (self.n * self.n)
        y = (pos_id % (self.n * self.n)) // self.n
        z = pos_id % self.n
        frac_x = x / self.n
        frac_y = y / self.n
        frac_z = z / self.n
        if frac:
            return frac_x, frac_y, frac_z,self.inverse_id(atom_id)
        else :
            return x, y, z, self.inverse_id(atom_id)

    # GRABBER FUNCTIONS
    def grab_positions(self, atom_id):
        """
        :param type:
        :return array[]

        - Allows us to grab all the int IDs associated with a given type
        - These essentially are just grabbing all positions
        """
        return[
            self.encode_var(x,y,z,atom_id, pos_rounding="int",frac = False)
            for x in range(self.n)
            for y in range(self.n)
            for z in range(self.n)
        ]

    def grab_types(self, x,y,z, pos_rounding, frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param pos_rounding:
        :param frac:
        :return: array[] of ints

        - grabs all variable IDs associated with a given position
        - can use fractional or integer positions
        - rounding incredibly useful for non valid positions!
        """
        return [self.encode_var(x,y,z,t, pos_rounding = pos_rounding ,frac = frac) for t in range(self.k)]

    def grab_positions_by_symbol(self, symbol):
        return self.grab_positions(self.atom_id(symbol))

    # NEIGHBOUR FUNCTIONS
    def create_offsets(self, distance):
        """
        :param distance:
        :return: tuple of distances

        - Calculates offsets based on distance
        - Main helper function for create_neighbours
        """
        offsets = []

        for dx,dy,dz in product(range( -distance, distance+1), repeat = 3):
            distance_calculated =  sqrt(dx**2 + dy**2 + dz**2)
            if distance_calculated == 0 :
                continue
            if distance_calculated == distance:
                 offsets.append((dx,dy,dz))
        return offsets


    def create_neighbours(self, x,y,z,distance, pos_rounding, dis_rounding, frac ):
        """
        :param x:
        :param y:
        :param z:
        :param distance:
        :param pos_rounding:
        :param dis_rounding:
        :param frac:
        :return: array[] of tuples nx,ny,nz

        - Allows us to create neighbours for a given position and given distance
        - pos_rounding allows us to use fractional positions
        - dis_rounding allows us to use fractional distances

        LIMITS:
        - current model limited to only using fractional or integer systems
        - CANNOT use integer positions and fractional distances or vice verse
        """

        x_unit = y_unit = z_unit = None
        scaled_distance = None
        neighbours = []

        if frac:
            x_unit, y_unit, z_unit = self.get_unit(x,y,z,pos_rounding = pos_rounding)
            if dis_rounding== "floor":
                scaled_distance = int(distance * self.n)
            elif dis_rounding == "ceil":
                scaled_distance = ceil(distance * self.n)
            elif dis_rounding == "round":
                scaled_distance = int(round(distance*self.n))
        elif not frac and pos_rounding == "int" and dis_rounding =="int":
            scaled_distance = distance
            x_unit = x
            y_unit = y
            z_unit = z


        for dx,dy,dz in self.create_offsets(scaled_distance):

            nx, ny, nz = x_unit + dx, y_unit + dy, z_unit + dz
            if 0 <= nx < self.n and 0 <= ny < self.n and 0 <= nz < self.n:
                neighbours.append((nx,ny,nz))
        return neighbours


    def encode_neighbors(self, x,y,z,atom_id,distance,pos_rounding,dis_rounding, frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param atom_id:
        :param distance:
        :param pos_rounding:
        :param dis_rounding:
        :param frac:
        :return: array[] of 6 tuple nx,ny,nz,type,pos_rounding,frac

        - returns neighbours based on a rounded distance along with a given type
        - dis_rounding is strictly used to determine nx,ny,nz
        - the 6 tuple can then be passed onto encode_var to find the ID
        - useful for creating constraints

        EXAMPLE:

        neighbors = encode_neighbors(0,34,0,63,0.79,symbol_id("Fe"),0.16,"floor","round",True)
        gives us an array of 6 tuples which corresponds:

        - all neighbours to x,y,z = floored (0.34,0.63,0.79)
        - which have a distance d = rounded(0.16)
        - all have type "Fe"
        - all int form perfect for encode_vars

        """
        neighbors_vars = []
        for nx,ny,nz in self.create_neighbours(x,y,z, distance,pos_rounding,dis_rounding, frac = frac ):

            neighbors_vars.append((nx,ny,nz,atom_id,pos_rounding,frac))

        return neighbors_vars

    # SETUP FUNCTION


    def initialise(self):
        """
        :return: constraints
        - Initialises all the variables in the crystal
        - forces so that at most one atom can be at each *valid* x,y,z position
        - only works with ints as these *always* corrsepond to *valid* sites

        """
        for x in range(self.n):
            for y in range(self.n):
                for z in range(self.n):
                    types_for_this_pos =self.grab_types(x,y,z,pos_rounding="int",frac = False)
                    for i in range(self.k):
                        for j in range(i+1,self.k):
                            self.cnf.append([-types_for_this_pos[i],-types_for_this_pos[j]])

    # FORCE CONSTRAINT FUNCTIONS

    def force_atom_at_position(self, x,y,z,atom_id,pos_rounding,frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param atom_id:
        :param pos_rounding:
        :param frac:
        :return: CONSTRAINTS

        - Allows user to force an atom of a certain type to a given position
        - Works for both fractional and integer systems
        """
        self.cnf.append([self.encode_var(x,y,z,atom_id,pos_rounding = pos_rounding,frac = frac)])

    def forbid_atom_at_position(self, x,y,z,atom_id,pos_rounding,frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param atom_id:
        :param pos_rounding:
        :param frac:
        :return: CONSTRAINTS

        - Allows user to forbid an atom of a certain type to a given position
        - Works for both fractional and integer systems
        """
        self.cnf.append([-self.encode_var(x,y,z,atom_id,pos_rounding = pos_rounding,frac = frac)])

    def require_one_of_types_at_position(self, x,y,z,atom_ids,pos_rounding ,frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param atom_ids: array []
        :param pos_rounding:
        :param frac:
        :return: CONSTRAINTS
        - Forces one of chosen type to be at a certain position
        - works with fractional and integer systems
        """
        vars = [self.encode_var(x,y,z,atom_ids,pos_rounding = pos_rounding, frac = frac) for atom_id in atom_ids]
        self.cnf.append(vars)

    def forbid_types_at_position(self, x, y, z, atom_ids,pos_rounding,frac = True):
        """
        :param x:
        :param y:
        :param z:
        :param atom_ids:
        :param pos_rounding:
        :param frac:
        :return: CONSTRAINTS
        - Builds on top of forbid_atom_at_position
        - can forbit multiple types at specific position
        - accepts fractional and integer systems
        """
        for atom_id in atom_ids:
            self.forbid_atom_at_position(x,y,z,atom_id, pos_rounding = pos_rounding, frac = frac)


    # ISOLATION FUNCTIONS

    def isolate_from_types(self, target_id, forbidden_neighbor_ids, distance, dis_rounding, frac = True):
        """
        :param target_id:
        :param forbidden_neighbor_ids: array []
        :param distance:
        :param dis_rounding:
        :param frac:
        :return: CONSTRAINTS
        - Allows user to isolate a certain type at all positions within a distance
        - from an array[] of forbidden types
        - only uses dis_rounding as it isolates from *valid* sites
        """
        for x in range(self.n):
            for y in range(self.n):
                for z in range(self.n):

                    atom = self.encode_var(x,y,z,target_id,pos_rounding ="int",frac = False)

                    for neighbor_id in forbidden_neighbor_ids:
                        if frac:
                            frac_x , frac_y , frac_z = self.get_frac(x,y,z)

                            neighbors = self.encode_neighbors(
                            frac_x, frac_y, frac_z,
                            neighbor_id,
                            distance = distance,
                            pos_rounding = "floor", #encodes valid site so does not matter
                            dis_rounding = dis_rounding,
                            frac = frac
                            )
                            for neighbor in neighbors:

                                self.cnf.append([-atom,-self.encode_var(*neighbor)])

                        if not frac and dis_rounding == "int":
                            # this is wrong?
                            # this pulls all nei
                            neighbors = self.encode_neighbors(x,y,z,neighbor_id, distance ,pos_rounding = "int",dis_rounding ="int", frac = False)
                            for neighbor in neighbors:
                                self.cnf.append([-atom, -self.encode_var(*neighbor)])


    def isolate(self,atom_id, distance,dis_rounding, frac = True):
        self.isolate_from_types(atom_id,list(range(self.k)),distance, dis_rounding = dis_rounding, frac = frac)

    def isolate_from_itself(self, atom_id, distance, dis_rounding, frac = True):
        forbidden_neigbour_types = [atom_id]
        self.isolate_from_types(atom_id,forbidden_neigbour_types,distance, dis_rounding = dis_rounding, frac = frac)

    # BY SYMBOL functions
    def isolate_from_types_by_symbol(self, target_symbol, forbidden_symbols, distance, dis_rounding, frac = True):
        target_type = self.atom_id(target_symbol)
        forbidden_neighbor_types = [self.atom_id(sym) for sym in forbidden_symbols]
        self.isolate_from_types(target_type, forbidden_neighbor_types, distance, dis_rounding = dis_rounding, frac = frac)

    def isolate_by_symbol(self, symbol, distance, dis_rounding, frac = True):
        self.isolate(self.atom_id(symbol), distance, dis_rounding = dis_rounding, frac = frac)

    def isolate_from_itself_by_symbol(self, symbol,distance, dis_rounding, frac = True):
        self.isolate_from_itself(self.atom_id(symbol),distance, dis_rounding = dis_rounding, frac = frac)



    # GLOBAL CONSTRAINTS :

    def get_forced(self,atom_id):

        """
        :param atom_id:
        :return: int

        returns all the vars that have been forced true for a specific type
        helper function for global cardinality constraints
        """
        forced_true_vars = [
            clause[0] for clause in self.cnf.clauses
            if len(clause) == 1 and clause[0] > 0 and (clause[0] - 1) % self.k == atom_id

        ]
        return forced_true_vars


    def at_most(self, amount, atom_id):

        # grabs all vars of atom_ type type
        vars = self.grab_positions(atom_id)
        atom_type = self.inverse_id(atom_id)

        forced_vars = self.get_forced(atom_id = atom_id)

        # adjust limit for SAT
        forced_count = len(forced_vars)
        adjusted_limit = amount - forced_count
        if adjusted_limit < 0:
            raise ValueError(
                f"Already forced {forced_count} atoms of type {atom_type} with atom ID {atom_id}"
                f" which exceeds at_most({amount})"
            )

        #Apply cardinality constraint to unforced variables only
        unforced_vars = [v for v in vars if v not in forced_vars ]
        if unforced_vars and adjusted_limit > 0:
            card = CardEnc.atmost(lits=unforced_vars, bound=adjusted_limit, encoding=EncType.totalizer)
            self.cnf.extend(card.clauses)


    def at_least(self, amount, atom_id):
        vars = self.grab_positions(atom_id)
        atom_type = self.inverse_id(atom_id)
        forced_vars = self.get_forced(atom_id = atom_id)
        forced_count = len(forced_vars)
        adjusted_bound = amount - forced_count
        if adjusted_bound <= 0:
            print(
                f"Note: Already forced {forced_count} atoms of type {atom_type} with atom ID {atom_id},"
                f"satisfying at_least({amount})"
            )
            return

        unforced_vars = [v for v in vars if v not in forced_vars ]
        if unforced_vars and adjusted_bound > 0 :
            card = CardEnc.atleast(lits = unforced_vars,
            bound= adjusted_bound, encoding=EncType.totalizer)
            self.cnf.extend(card.clauses)

    def strictly(self, amount, atom_id):
        all_vars = self.grab_positions(atom_id)

        # Detect already forced True variables
        forced_vars = self.get_forced(atom_id = atom_id)

        forced_count = len(forced_vars)

        adjusted_lower = amount - forced_count
        adjusted_upper = amount - forced_count

        if adjusted_lower < 0 or adjusted_upper < 0:
            raise ValueError(
                f"Already forced {forced_count} atoms of type {atom_id}, which violates strictly({amount})."
            )

        # Apply constraints to unforced variables
        unforced_vars = [v for v in all_vars if v not in forced_vars]
        if unforced_vars:
            if adjusted_lower > 0:
                card1 = CardEnc.atleast(lits=unforced_vars, bound=adjusted_lower, encoding=EncType.totalizer)
                self.cnf.extend(card1.clauses)
            if adjusted_upper >= 0:
                card2 = CardEnc.atmost(lits=unforced_vars, bound=adjusted_upper, encoding=EncType.totalizer)
                self.cnf.extend(card2.clauses)

    # Solves model
    def solve(self, solver_name="glucose3"):
        from pysat.solvers import Solver

        with Solver(name=solver_name, bootstrap_with=self.cnf.clauses) as solver:
            is_sat = solver.solve()
            if is_sat:
                return solver.get_model()
            else:
                return None

    # Decodes solution into readable format
    # does not decode correctly
    def decode_solution(self ,model, frac = True):
        max_original = self.n ** 3 * self.k
        variables = []
        for encoded_var in model:
           # truth
            if 0<encoded_var<= max_original:
                x,y,z, atom_symbol = self.inverse_var(encoded_var,frac = frac)
                variables.append((x,y,z,atom_symbol,True))
        return variables


    # export to CIF

    def export_to_cif(self, model, filename="output.cif"):
        """
        Export a solve SAT model to a CIF file
        """
        #Defines a cubic lattice
        lattice = Lattice.cubic(self.n)

        #Extract atoms
        species = []
        frac_coordinates = []
        decoded = self.decode_solution(model)
        for x,y,z,symbol,truth_value in decoded:
            if truth_value:
                fx = x/self.n
                fy = y/self.n
                fz = z/self.n
                species.append(symbol)
                frac_coordinates.append([fx,fy,fz])

        structure = Structure(lattice,species,frac_coordinates)
        structure.to(filename)
        print(f"CIF file is saved to {filename}")








