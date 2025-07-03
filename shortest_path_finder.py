import networkx as nx
import itertools
from rhg import *


class ShortestPathFinder:
    def __init__(self, G):
        self.G = G

    def dijkstra(self, start, end):
        """
        Find the shortest path between two nodes in the graph represented by adjacency matrix A.
        :param start: Starting node index.
        :param end: Ending node index.
        :return: List of nodes in the shortest path, or None if no path exists.
        """
        try:
            path = nx.shortest_path(self.G, source=start, target=end)
            return path
        except nx.NetworkXNoPath:
            print("No path exists between the specified nodes.")
            return None
    
    def hyperbolic_geodesic_distance(self, start, end, rtheta_coords, length):
        middle_length = length - 2
        geodesic_distance = dict()
        for i in self.G.nodes:
            if i != start and i != end:
                d_is = rhg.distance(rtheta_coords[i], rtheta_coords[start])
                d_it = rhg.distance(rtheta_coords[i], rtheta_coords[end])
                d_st = rhg.distance(rtheta_coords[start], rtheta_coords[end])
                geodesic_distance[i] = 1/2 * (d_is + d_it - d_st) + np.log(2)
        sorted_distances = sorted(geodesic_distance.items(), key=lambda x: x[1])  # sort by distance
        candidate_nodes_coords = {node:rtheta_coords[node] for node, _ in sorted_distances[:middle_length]}  # get the top 'length' nodes
        sorted_path = sorted(candidate_nodes_coords.items(), key=lambda x: rhg.distance(x[1], rtheta_coords[start]))  # sorted by geodesic distance to start
        return [start] + [node for node, _ in sorted_path] + [end]
    
    def dijkstra_weighted_rhg_p(self, start, end, rtheta_coords, R, T):
        nodes = self.G.nodes
        for i, j in itertools.combinations(nodes, 2):
            if not self.G.has_edge(i, j) and not self.G.has_edge(j, i):
                delta_theta = np.pi-np.abs(np.pi-np.abs(rtheta_coords[i][1]-rtheta_coords[j][1]))
                d = np.arccosh(np.cosh(rtheta_coords[i][0])*np.cosh(rtheta_coords[j][0])-np.sinh(rtheta_coords[i][0])*np.sinh(rtheta_coords[j][0])*np.cos(delta_theta))
                p = 1/(np.exp((d-R)/T)+1)
                self.G.add_edge(i, j, weight=1/p)
        path = nx.shortest_path(self.G, source=start, target=end, weight="weight")
        return path

    
    @staticmethod
    def precision(predicted_path, actual_paths):
        """
        Calculate the precision of the predicted path against the actual path.
        :param predicted_path: List of nodes in the predicted path.
        :param actual_paths: List of List of nodes in the actual path.
        :return: Precision value as a float.
        """
        if not predicted_path or not actual_paths:
            return 0.0
        true_positives = 0
        actual_shortest_path_length = None
        for actual_path in actual_paths:
            if actual_shortest_path_length is None:
                actual_shortest_path_length = len(actual_path) - 2  # Exclude start and end nodes
            joint_middle_nodes_number = len(set(predicted_path[1:-1]) & set(actual_path[1:-1]))
            true_positives = joint_middle_nodes_number if joint_middle_nodes_number > true_positives else true_positives
        return true_positives / actual_shortest_path_length
    
    @staticmethod
    def overlap(predicted_path, actual_paths):
        """
        Calculate the overlap of the predicted path against the actual paths.
        :param predicted_path: List of nodes in the predicted path.
        :param actual_paths: List of List of nodes in the actual paths.
        :return: Overlap value as a float.
        """
        if not predicted_path or not actual_paths:
            return -1
        max_jackard_index = 0
        for actual_path in actual_paths:
            joint_middle_nodes_number = len(set(predicted_path[1:-1]) & set(actual_path[1:-1]))
            if len(set(predicted_path[1:-1])) == 0 and len(set(actual_path[1:-1])) == 0:
                jackard_index = 1
            else:
                jackard_index = joint_middle_nodes_number / len(set(predicted_path[1:-1]) | set(actual_path[1:-1]))
            max_jackard_index = max(max_jackard_index, jackard_index)
        return max_jackard_index
    
    @staticmethod
    def path_edit_distance(predicted_path, actual_path):
        """
        Calculate the edit distance between two paths using dynamic programming.
        
        The edit distance is the minimum number of operations (insertions, deletions, 
        or substitutions) needed to transform one path into another.
        
        Args:
            predicted_path: List of nodes representing the first path
            actual_path: List of nodes representing the second path
            verbose: If True, prints the DP table and operations
        
        Returns:
            int: The minimum edit distance between the two paths
        """
        m, n = len(predicted_path), len(actual_path)
        # Create DP table - dp[i][j] represents edit distance between
        # predicted_path[:i] and actual_path[:j]
        dp = [[0] * (n + 1) for _ in range(m + 1)]
        # Initialize base cases
        # Converting empty path to actual_path[:j] requires j insertions
        for j in range(n + 1):
            dp[0][j] = j
        # Converting predicted_path[:i] to empty path requires i deletions
        for i in range(m + 1):
            dp[i][0] = i
        # Fill the DP table
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if predicted_path[i-1] == actual_path[j-1]:
                    # No operation needed if nodes are the same
                    dp[i][j] = dp[i-1][j-1]
                else:
                    # Take minimum of three operations:
                    # 1. Insert: dp[i][j-1] + 1
                    # 2. Delete: dp[i-1][j] + 1  
                    # 3. Substitute: dp[i-1][j-1] + 1
                    dp[i][j] = min(
                        dp[i][j-1] + 1,    # insertion
                        dp[i-1][j] + 1,    # deletion
                        dp[i-1][j-1] + 1   # substitution
                    )
        return dp[m][n]
    
    @staticmethod
    def path_edit_similarity(predicted_path, actual_paths):
        if not predicted_path or not actual_paths:
            return -1
        max_similarity = 0
        for actual_path in actual_paths:
            edit_dist = ShortestPathFinder.path_edit_distance(predicted_path, actual_path)
            max_len = max(len(predicted_path), len(actual_path))
            if max_len == 0:
                similarity =  1  # Both paths are empty
            else:
                similarity = 1 - (edit_dist / max_len)
            max_similarity = max(max_similarity, similarity)
        return max_similarity