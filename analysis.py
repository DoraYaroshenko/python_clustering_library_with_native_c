from sklearn.metrics import silhouette_score
import sys
import numpy as np
import pandas as pd
import math
import symnmf
EPS = 0.0001
ITER = 300

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

# oboze 
def kmeans(K,vectors):
    points = vectors.values.tolist()
    centroids = []
    for i in range(K):
        centroids.append(points[i].copy())
    for i in range(ITER):
        assigned_to_cent = [ [] for _ in range(len(centroids)) ]
        convergence = True
        for i in points: # i in both loops??? 
            min_dist_from_centroid = sys.float_info.max
            index_of_centroid_with_min_dist = 0
            for j in range(K):
                centroid = centroids[j]
                distance = euclidian_distance(i, centroid)
                if (distance<min_dist_from_centroid):
                    min_dist_from_centroid = distance
                    index_of_centroid_with_min_dist = j
            assigned_to_cent[index_of_centroid_with_min_dist].append(i)
    
        for j in range(len(centroids)):
            if not assigned_to_cent[j]:
                continue
            mean = [0.0 for x in range(len(centroids[j]))]
            for assigned_centroid in assigned_to_cent[j]:
                for i in range(len(assigned_centroid)):
                    mean[i] += assigned_centroid[i]
            for i in range(len(mean)):
                mean[i] /= len(assigned_to_cent[j])
            if (euclidian_distance(centroids[j], mean))>=EPS:
                convergence = False
            centroids[j] = mean
        if convergence:
            break
    return centroids

def getInputVariables():
    if len(sys.argv) != 3:
        # errorHandling! 
        print("An Error Has Occurred")
        sys.exit(1)
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
    K, file_name = getInputVariables()
    vectors = symnmf.getDataPoints(file_name)
    N = len(vectors)
    if (K>=N):
        print("Incorrect number of clusters!")
        sys.exit(1)
    points_array=vectors.to_numpy()
    H = symnmf.getFinalH(points_array,K)
    labels_sym = np.argmax(H, axis=1)
    sym_score = silhouette_score(points_array, labels_sym)
    final_centroids = kmeans(K,vectors)
    labels_kmeans = [get_cluster_association(vec, final_centroids) for vec in points_array]
    kmeans_score = silhouette_score(points_array, labels_kmeans)
    print(f"nmf: {sym_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")