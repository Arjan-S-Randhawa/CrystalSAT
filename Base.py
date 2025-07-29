
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

        # specifying geometry of unit cell
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.allowed = allowed
        self.cnf = CNF()
        self.allowed = allowed
        self.use_allowed = bool(allowed)
        self.k = len(allowed) if self.use_allowed else 118
        self.max_real = self.n_x * n_y * n_z * self.k
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
        self.populate_var_dict()
