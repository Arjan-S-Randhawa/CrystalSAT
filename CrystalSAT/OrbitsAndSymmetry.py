from pymatgen.symmetry.groups import SpaceGroup
import numpy as np
from pysat.card import CardEnc

class OrbitsAndSymmetryMixin:

    def _canon_coord_triplet(self, arr, eps=1e-9, ndp=8):
        """Return a (x,y,z) tuple of plain Python numbers suitable for hashing.
        - Wrap to [0,1)
        - Round to ndp decimals
        - Convert near-integers to int for pretty tuples like (0, 0, 0)
        """
        a = np.mod(np.asarray(arr, dtype=float), 1.0)
        a = np.round(a, ndp)
        out = []
        for v in a:
            rv = round(float(v))
            if abs(float(v) - rv) < eps:
                out.append(int(rv) % 1)  # 1 -> 0 after mod
            else:
                out.append(float(v))
        return tuple(out)


    def get_orbit_positions(self, x, y, z, space_group, tol=1e-5):
        """
        Returns the symmetry positions for a given x, y, z in a space group.
        :param x: x coordinate
        :param y: y coordinate
        :param z: z coordinate
        :param space_group: space group number or symbol
        :return: list of symmetry positions
        """
        sg = SpaceGroup(space_group)

        # Create the input position as fractional coordinates
        position = np.array([x, y, z])

        # Get orbit (all symmetry equivalent positions)
        orbit = sg.get_orbit(position, tol=tol)
        symmetry_positions = []
        for pos in orbit:
            symmetry_positions.append(self._canon_coord_triplet(pos))
        return symmetry_positions


    def populate_orbit_dict(self, space_group, tol=1e-5):
        """
        Populates the orbit dictionary with symmetry positions for each (x, y, z).
        :param space_group: space group number or symbol
        :param tol: tolerance for symmetry operations
        :return: None
        """
        self.orbits = []
        seen = set()
        if space_group is None:
            return # No space group provided, no orbits to populate
        else:

            for x, y, z in self.positions:
                x_frac,y_frac,z_frac = self.to_frac(x, y, z, system="int", pos_rounding="round")
                orbit_positions = self.get_orbit_positions(x_frac, y_frac, z_frac, space_group, tol)
                canonical = tuple(sorted(orbit_positions))
                if canonical in seen:
                    continue
                seen.add(canonical)
                orbit_id = len(self.orbits)
                self.orbits.append(list(canonical))
                self.orbit_dict[orbit_id] = sorted(orbit_positions)

                for pos in sorted(orbit_positions):
                    x1,y1,z1 = pos
                    self.inverse_orbit_dict[(x1,y1,z1)] = orbit_id


    def orbit_var(self, orbit_id, atom_id):
        if self.space_group is None:
            raise ValueError("No space group defined, cannot create orbit variable.")

        return self.orbit_pool.id(("orbit", orbit_id, atom_id))

    def choose_orbits(self, min_count, max_count, atom_id, eligible_orbits=None, weight=None):
        """
        Constrain the solver to choose exactly one orbit for this atom_id
        and link that orbit variable to all its site variables.
        If eligible_orbits is provided, restrict the choice to that subset.
        """

        forced_orbits_count = len(self.grab_forced_orbits(atom_id))
        if eligible_orbits is None:
            orbit_vars = [self.orbit_var(oid, atom_id) for oid in self.orbit_dict.keys()]
            if min_count is not None:
                adjusted_min_count = min_count - forced_orbits_count
                if adjusted_min_count < 0:
                    raise ValueError(f"min_count {min_count} is less than the number of forced orbits {forced_orbits_count} for atom type {atom_id}.")
                elif adjusted_min_count > forced_orbits_count + len(self.orbit_dict):
                    raise ValueError(f"min_count {min_count} is greater than the number of available orbits {len(self.orbit_dict)} + forced orbits {forced_orbits_count} for atom type {atom_id}.")

                atleast = CardEnc.atleast(lits=orbit_vars, bound=adjusted_min_count, encoding=1, vpool=self.vpool)

                for cl in atleast.clauses:
                    if weight is None:
                        self.cnf.append(cl)
                    else:
                        self.cnf.append(cl, weight=weight)

                self.vpool.start_from = self.vpool.id

            if max_count is not None:
                adjusted_max_count = max_count - forced_orbits_count
                if adjusted_max_count < 0:
                    raise ValueError(f"max_count {max_count} is less than the number of forced orbits {forced_orbits_count} for atom type {atom_id}.")
                elif adjusted_max_count > len(self.orbit_dict):
                    raise ValueError(f"max_count {max_count} is greater than the number of available orbits {len(self.orbit_dict)} + forced orbits {forced_orbits_count} for atom type {atom_id}.")
                atmost = CardEnc.atmost(lits=orbit_vars, bound=adjusted_max_count, encoding=1, vpool=self.vpool)
                for cl in atmost.clauses:
                    if weight is None:
                        self.cnf.append(cl)
                    else:
                        self.cnf.append(cl, weight=weight)

                self.vpool.start_from = self.vpool.id



        for orbit in self.orbit_dict.keys():
            orbit_var = self.orbit_var(orbit, atom_id)
            positions = self.orbit_dict[orbit]
            for x, y, z in positions:
                # Link the orbit variable to all its sites
                x_int, y_int, z_int = self.to_int(x,y,z,system="frac", pos_rounding="round")
                site_var = self.encode_var(x_int, y_int, z_int, atom_id)
                self.cnf.append([-orbit_var, site_var])







    # turn to weighted later
    def force_orbit(self, atom_id, orbit_id):
        """
        Force a specific orbit for an atom_id.
        This will create a variable for the orbit and link it to all its site variables.
        :param atom_id: ID of the atom type
        :param orbit_id: ID of the orbit to force
        """
        if self.space_group is None:
            raise ValueError("No space group defined, cannot force orbit.")

        if orbit_id not in self.orbit_dict:
            raise ValueError(f"Orbit ID {orbit_id} not found. Did you call populate_orbit_dict()?")

        # Get SAT var for this (orbit, atom) and assert it true
        orbit_var = self.orbit_var(orbit_id, atom_id)
        self.cnf.append([orbit_var])  # force the orbit to be active

        # Link orbit var to all symmetry-equivalent sites (and conversely)
        for x, y, z in self.orbit_dict[orbit_id]:

            x_int, y_int, z_int = self.to_int(x,y,z,system="frac", pos_rounding="round")

            self.force_atom_at_position(x_int, y_int, z_int, atom_id)


    def grab_forced_orbits(self,atom_id):

        forced = []

        if hasattr(self.cnf, "hard") and hasattr(self.cnf, "soft"):
            # WCNF format
            all_clauses = self.cnf.hard
        else:
            # CNF format
            all_clauses = self.cnf.clauses

        for clause in all_clauses:
            if len(clause) == 1 and clause[0] > 0:
                var  = clause[0]
                obj = self.orbit_pool.obj(var)
                if obj and obj[0] == "orbit" and obj[2] == atom_id:
                    _,orbit_id,_ = obj
                    forced.append(orbit_id)

        return forced















