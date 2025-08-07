
from ase import Atoms
from ase.cell import Cell
from pysat.formula import CNF
from pysat.formula import IDPool

from .Encoding import EncodingMixin
from .Coordinate import CoordinateMixin
from .Get import GetMixin
from .NeighborAndDistances import NeighborAndDistancesMixin
from .Constraints import ConstraintsMixin
from .NeighborConstraints import NeighborConstraintsMixin
from .Grab import GrabMixin
from .Cardinality import CardinalityMixin
from .SolveAndExport import SolveAndExportMixin


class CrystalSAT(EncodingMixin,CoordinateMixin,GetMixin,
                 NeighborAndDistancesMixin,ConstraintsMixin,
                  NeighborConstraintsMixin,GrabMixin,CardinalityMixin,
                  SolveAndExportMixin):

    def __init__(self, n_x,n_y,n_z,
                       a,b,c,alpha,
                       beta,gamma,allowed):

        # specifying grid dimensions
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z

        # specifying unit cell parameters
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        self.allowed = allowed

        # changes lower bound and use_allowed accordingly
        if len(allowed) > 0:
            self.use_allowed = True
            self.lower = 0
        else:
            self.use_allowed = False
            self.lower = 1

        # dictionary of ions mapped to atom IDs
        self.ion_dict = {}
        self.reverse_ion_dict = {}
        self.populate_ion_dict()


        # k is the number of allowed atom types and now ions
        self.k =  len(allowed) if self.use_allowed else 118 + len(self.ion_dict)
        self.max_real = self.n_x * n_y * n_z * self.k

        # Initialize CNF and IDPool
        self.cnf = CNF()
        self.vpool = IDPool(start_from=self.max_real + 1)

        # ASE cell object
        self.cell = Cell.fromcellpar([self.a,self.b,self.c,self.alpha,self.beta,self.gamma])
        self.valid_positions = [ (x,y,z) for x in range(self.n_x) for y in range(self.n_y) for z in range(self.n_z) ]
        symbols = ['X'] * len(self.valid_positions)
        frac_positions = [ (x/self.n_x, y/self.n_y, z/self.n_z) for (x,y,z) in self.valid_positions ]
        cart_positions = [ self.cell.cartesian_positions([fp])[0] for fp in frac_positions ]

        # Creates grid used for geometry calculations
        self.grid = Atoms(symbols=symbols, positions=cart_positions, cell=self.cell, pbc=True)

        # Populate variable dictionary after all attributes are set
        self.var_dict = {}
        self.reverse_var_dict = {}
        self.populate_var_dict()

        # List of all positions in the grid
        self.positions = [(x, y, z) for x in range(self.n_x) for y in range(self.n_y) for z in range(self.n_z)]
