from sklearn.metrics import silhouette_score
import sys
import numpy as np
import pandas as pd
import math
import symnmf
EPS = 0.0001
MAX_ITERATIONS = 300

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

def assign_to_closest_cluster(points,centroids,K):
    clusters_with_assigned_points = [ [] for _ in range(K) ]
    for point in points:
        min_dist_from_centroid = sys.float_info.max
        index_of_centroid_with_min_dist = 0
        for i in range(K):
            centroid = centroids[i]
            distance = euclidian_distance(point, centroid)
            if (distance<min_dist_from_centroid):
                min_dist_from_centroid = distance
                index_of_centroid_with_min_dist = i
        clusters_with_assigned_points[index_of_centroid_with_min_dist].append(point)
    return clusters_with_assigned_points

def updateCentroids(centroids, clusters_with_assigned_points, K, dimension):
    convergence = True
    for i in range(K):
        if not clusters_with_assigned_points[i]:
            continue
        num_of_assigned_points = len(clusters_with_assigned_points[i])
        new_centroid = [0.0 for x in range(dimension)]
        for point in clusters_with_assigned_points[i]:
            for j in range(dimension):
                new_centroid[j] += point[j]
        for j in range(dimension):
            new_centroid[j] /= num_of_assigned_points
        if (euclidian_distance(centroids[i], new_centroid))>=EPS:
            convergence = False
        centroids[i] = new_centroid
    return centroids,convergence

def kmeans(K,vectors):
    points = vectors.values.tolist()
    dimension = len(points[0])
    centroids = []
    for i in range(K):
        centroids.append(points[i].copy())
    for i in range(MAX_ITERATIONS):
        clusters_with_assigned_points = assign_to_closest_cluster(points,centroids,K)
        centroids,convergence = updateCentroids(centroids,clusters_with_assigned_points,K,dimension)
        if convergence:
            break
    return centroids

def get_input_variables():
    if len(sys.argv) != 3:
        symnmf.error_handling()
    K = int(sys.argv[1])
    file_name = sys.argv[2]
    return K, file_name

def get_cluster_association(vector, centroids):
    min_delta = float("inf")
    association = 0
    for i, centroid in enumerate(centroids):
        delta = euclidian_distance(vector, centroid)
        if delta < min_delta:
            min_delta = delta
            association = i
    return association

if __name__ == "__main__":
    K, file_name = get_input_variables()
    vectors = symnmf.get_data_points(file_name)
    N = len(vectors)
    if (K>=N):
        print("Incorrect number of clusters!")
        sys.exit(1)
    points_array=vectors.to_numpy()
    H = symnmf.calculate_decomposition_matrix(points_array,K)
    labels_sym = np.argmax(H, axis=1)
    sym_score = silhouette_score(points_array, labels_sym)
    final_centroids = kmeans(K,vectors)
    labels_kmeans = [get_cluster_association(vec, final_centroids) for vec in points_array]
    kmeans_score = silhouette_score(points_array, labels_kmeans)
    print(f"nmf: {sym_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")