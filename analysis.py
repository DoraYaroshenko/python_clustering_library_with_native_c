from sklearn.metrics import silhouette_score
import sys
import numpy as np
import symnmfmodule
import mykmeanssp
import pandas as pd
import math

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

def similarity_matrix(points_array):
    squares_sum = np.sum(points_array**2, axis=1, keepdims=True) # לכל שורה במטריצה מעלה את האיבר בכל עמודה בריבוע וסוכם את העמודות. זה נהיה וקטור עמודה שמכיל את הסכומים הנ״ל בכל שורה
    squared_sub_mat = squares_sum+squares_sum.T # זה מטריצה שבה xij = |x_i|^2+|x_j|^2
    squared_sub_mat = squared_sub_mat-2*(points_array @ points_array.T) # מכל איבר איקס איי ג׳יי במטריצה הקודמת יצא לנו להחסיר שתיים כפול מ״פ בין איקס איי לאיקס ג׳יי 
    # סהכ מתקיים מ״חוק הקוסינוסים״ ש*לא* למדנו עם סמיון: ||x_i - x_j||^2 = ||x_i|| + ||x_j|| - 2*<x_i, x_j> and it matches what we have done so far
    A = np.exp(-squared_sub_mat/2) # מטריצה לפי ההגדרה, עד כדי זה שהאלכסון לא תואם להגדרה
    np.fill_diagonal(A, 0) # איפוס האלכסון
    return A

def degree_matrix(A):
    d = np.sum(A, axis=1) # סוכם את האיברים בכל שורה
    D = np.diag(d) # מציב את האיברים של די, שזה הדרגות שחישבנו קודם, על האלכסון - ושאר המטריצה היא אפסים
    return D

def normalized_similarity_mat(D, A):
    d = np.diag(D)
    D_power = np.diag(1.0 / np.sqrt(d)) # מעלים את האלכסון בחזקת מינוס חצי שזה כמו אחד חלקי שורש ומציבים אותו באלכסון המטריצה החדשה שלנו
    W = D_power @ A @ D_power # מחשבים את הביטוי שהוגדר
    return W

def init_H(W, K):
    n = W.shape[0]
    m = W.mean()
    upper = 2*np.sqrt(m/K)
    H = np.random.uniform(0, upper, size=(n, K))
    return H

def getInputVariables():
    if len(sys.argv) != 3:
        print("An Error Has Occurred")
        sys.exit(1)
    K = int(sys.argv[1])
    file_name = sys.argv[2]
    return K, file_name

def getObservations(file_name):
    vectors = pd.read_csv(file_name, header=None)
    columns_names = [f"coordinate {i}" for i in range(vectors.shape[1])]
    vectors.columns = columns_names
    return vectors

def initCentroids(vectors,points_array,K):
    points = np.copy(points_array)
    np.random.seed(1234)
    index_chosen = np.random.choice(len(points_array))
    centroid = points[index_chosen]
    centroids = np.zeros((K,len(points[0])))
    centroids_indexes = []
    centroids[0]=centroid
    centroids_indexes.append(vectors.index[index_chosen])
    points = np.delete(points,index_chosen,axis=0)
    # print(points)
    for i in range(1,K):
        distances = []
        for p in points:
            # print(p)
            # print(centroids)
            point_distances = [euclidian_distance(centroid,p) for centroid in centroids[:i]]
            dist = min(point_distances)
            distances.append(dist)
        probabilities=[distance/sum(distances) for distance in distances]
        index_chosen=np.random.choice(len(points),p=probabilities)
        centroids[i]=points[index_chosen]
        centroids_indexes.append(vectors.index[np.where((points_array==centroids[i]).all(axis=1))[0][0]])
        points = np.delete(points,index_chosen, axis=0)
    return centroids,centroids_indexes

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
    vectors = getObservations(file_name) #ok
    N = len(vectors) #ok
    points_array=vectors.to_numpy() # okay, returns actual matrix in numpy form..

    A = similarity_matrix(points_array) # maybe should use sym method from c file instead.
    D = degree_matrix(A) # maybe should use ddg method from c file instead
    W = normalized_similarity_mat(D, A) # should have used norm from c file instead apparently
    H = init_H(W, K)
    H = symnmfmodule.symnmf(H, W)
    
    labels_sym = np.argmax(H, axis=1)
    sym_score = silhouette_score(points_array, labels_sym)

    centroids,centroids_indexes = initCentroids(vectors,points_array, K)
    final_centroids = mykmeanssp.fit(centroids, points_array, 300, 1e-4)
    labels_kmeans = [get_cluster_association(vec, final_centroids) for vec in points_array]
    kmeans_score = silhouette_score(points_array, labels_kmeans)

    print(f"KMeans Score: {kmeans_score:.4f}")
    print(f"Symnmf Score: {sym_score:.4f}")
