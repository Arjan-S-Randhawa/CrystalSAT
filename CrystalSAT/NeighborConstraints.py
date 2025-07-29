
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

        forbidden_neighbor_ids =[k for k in range(self.k) if k != target_id]
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

