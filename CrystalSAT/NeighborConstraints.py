
from itertools import combinations

class NeighborConstraintsMixin:


    def isolate_from_types(self, target_id, forbidden_neighbor_ids, cutoff, tolerance,ball= True):
        """
        Isolates a specific atom type from chosen types within a given distance.
        tolerance is used to allow for some flexibility in distance matching.
        :param target_id:
        :param forbidden_neighbor_ids: array []
        :param cutoff:
        :param tolerance:
        :param ball:
        :return:
        """
        for x in range(self.n_x):
            for y in range(self.n_y):
                for z in range(self.n_z):
                    atom = self.encode_var(x, y, z, target_id)
                    for neighbor_id in forbidden_neighbor_ids:
                        neighbors = self.encode_neighbors(x, y, z, atom_id= neighbor_id, cutoff = cutoff, tolerance = tolerance, system="int", pos_rounding="int", ball=ball)
                        for neighbor in neighbors:
                            self.cnf.append([-atom, -neighbor])


    def isolate(self, target_id, cutoff, tolerance, ball=True):
        """
        Isolates a specific atom type from all other types within a given distance.
        :param target_id:
        :param cutoff:
        :param tolerance:
        :param ball:
        :return:
        """

        forbidden_neighbor_ids =[k for k in range(self.lower,self.k) if k != target_id]
        self.isolate_from_types(target_id = target_id, forbidden_neighbor_ids = forbidden_neighbor_ids,cutoff= cutoff , tolerance=tolerance, ball=True)


    def isolate_from_itself(self, target_id, cutoff, tolerance, ball = True):
        """
        Isolates a specific atom type from itself within a given distance.
        :param target_id:
        :param cutoff:
        :param tolerance:
        :param ball:
        :return:
        """
        forbidden_neigbor_types = [target_id]
        self.isolate_from_types(target_id = target_id ,forbidden_neighbor_ids=forbidden_neigbor_types,cutoff = cutoff, tolerance = tolerance, ball=ball)


    def enforce_closest_dist(self, atom_id1, atom_id2, min_dist):

        """
        Enforces the closest distance between two atom types in the CNF.
        This means, atoms of type atom_id1 and atom_id2 must be at least min_dist apart.
        :param atom_id1:
        :param atom_id2:
        :param min_dist:
        :return:
        """
        for (x1, y1, z1), (x2, y2, z2) in combinations(self.positions, 2):

            idx1 = x1 * (self.n_y * self.n_z) + y1 * self.n_z + z1
            idx2 = x2 * (self.n_y * self.n_z) + y2 * self.n_z + z2

            rad_1 = self.get_radius(*self.inverse_id(atom_id1)) if isinstance(self.inverse_id(atom_id1),
                                                                              tuple) else self.get_radius(
                self.inverse_id(atom_id1))

            rad_2 = self.get_radius(*self.inverse_id(atom_id2)) if isinstance(self.inverse_id(atom_id2),
                                                                              tuple) else self.get_radius(
                self.inverse_id(atom_id2))

            dist = self.grid.get_distance(idx1, idx2, mic=True)

            if dist < rad_1 + rad_2 + min_dist:
                self.cnf.append([-self.encode_var(x1, y1, z1, atom_id1),
                                 -self.encode_var(x2, y2, z2, atom_id2)])

                self.cnf.append([-self.encode_var(x1, y1, z1, atom_id2),
                                 -self.encode_var(x2, y2, z2, atom_id1)])



