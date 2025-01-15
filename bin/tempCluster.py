#!/usr/bin/env python3

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

def validate_groups(distance_matrix, predefined_groups):
    within_distances = []
    between_distances = []
    n = len(predefined_groups)
    
    for i in range(n):
        for j in range(i + 1, n):
            if predefined_groups[i] == predefined_groups[j]:
                within_distances.append(distance_matrix[i, j])
            else:
                between_distances.append(distance_matrix[i, j])
    
    return np.mean(within_distances) < np.mean(between_distances)

# Example use
distance_matrix = np.array([[0, 1, 5], [1, 0, 4], [5, 4, 0]])  # example distance matrix
predefined_groups = [0, 0, 1]  # example group labels
validate_groups(distance_matrix, predefined_groups)







def assign_unknown_samples(distance_matrix, predefined_groups, threshold=1.0):
    n_samples = len(predefined_groups)
    new_groups = predefined_groups.copy()
    unknown_indices = [i for i in range(n_samples) if predefined_groups[i] == -1]  # -1 for unknown

    for idx in unknown_indices:
        distances = distance_matrix[idx]
        min_distance = np.min(distances)
        if min_distance < threshold:
            closest_group = predefined_groups[np.argmin(distances)]
            new_groups[idx] = closest_group
        else:
            new_groups[idx] = max(predefined_groups) + 1  # create a new group
    
    return new_groups

# Example use
predefined_groups = [0, 0, 1, -1, -1]  # with unknown samples
assign_unknown_samples(distance_matrix, predefined_groups)






def cluster_unknown_samples(distance_matrix, max_clusters=10, method='gmm'):
    if method == 'gmm':
        bic = []
        for n_components in range(1, max_clusters + 1):
            gmm = GaussianMixture(n_components=n_components)
            gmm.fit(distance_matrix)
            bic.append(gmm.bic(distance_matrix))
        best_n_components = np.argmin(bic) + 1
        gmm = GaussianMixture(n_components=best_n_components)
        gmm.fit(distance_matrix)
        return gmm.predict(distance_matrix)
    
    elif method == 'hierarchical':
        Z = linkage(squareform(distance_matrix), 'ward')
        clusters = AgglomerativeClustering(n_clusters=None, distance_threshold=0)
        return clusters.fit_predict(distance_matrix)

# Example use
cluster_unknown_samples(distance_matrix)
