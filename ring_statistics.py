import networkx as nx
import numpy as np

def construct_adjacency_matrix(num_atoms, coords, cutoff):
    adjacency_matrix = np.zeros((num_atoms, num_atoms))
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = np.linalg.norm(coords[i] - coords[j])
            if distance <= cutoff:
                adjacency_matrix[i, j] = 1
                adjacency_matrix[j, i] = 1
    return adjacency_matrix

def find_rings(num_atoms, coords, cutoff, max_ring_size):
    adjacency_matrix = construct_adjacency_matrix(num_atoms, coords, cutoff)
    G = nx.from_numpy_array(adjacency_matrix)
    all_rings = list(nx.cycle_basis(G))
    filtered_rings = [ring for ring in all_rings if len(ring) <= max_ring_size]
    return filtered_rings

def ring_statistics(num_atoms, coords, cutoff, max_ring_size):
    rings = find_rings(num_atoms, coords, cutoff, max_ring_size)
    return rings