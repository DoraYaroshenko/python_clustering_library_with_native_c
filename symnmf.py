import math
import sys
import numpy as np
import pandas as pd
import _symnmf

np.random.seed(1234)

def error_handling():
    print("An Error Has Occurred")
    sys.exit(1)

def init_H(W, K):
    W_as_array = np.array(W)
    n = W_as_array.shape[0]
    m = W_as_array.mean()
    upper = 2*np.sqrt(m/K)
    H = np.random.uniform(0, upper, size=(n, K))
    return H

def print_matrix(A):
    for row in A:
        line = ",".join(f"{val:.4f}" for val in row)
        print(line)

def get_input_variables():
    if len(sys.argv) != 4:
        error_handling()
    K = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return K, goal, file_name

def get_data_points(file_name):
    vectors = pd.read_csv(file_name, header=None)
    vectors.columns = [f"coordinate {i}" for i in range(vectors.shape[1])]
    vectors.index.name = "key"
    return vectors

def calculate_decomposition_matrix(points_array, K):
    W = _symnmf.norm(points_array)
    if W is None:
        error_handling()
    H = init_H(W, K)
    H = _symnmf.symnmf(H, W)
    if H is None:
        error_handling()
    return H

if __name__ == "__main__":
    K, goal, file_name = get_input_variables()
    vectors = get_data_points(file_name)
    N = len(vectors)
    if K>=N or K<=1:
        error_handling()
    points_array=vectors.to_numpy()
    match goal:
        case "sym":
            A = _symnmf.sym(points_array)
            if A is None:
                error_handling()
            print_matrix(A)
        case "symnmf":
            H = calculate_decomposition_matrix(points_array, K)
            print_matrix(H)
        case "ddg":
            D = _symnmf.ddg(points_array)
            if D is None:
                error_handling()
            print_matrix(D)
        case "norm":
            W = _symnmf.norm(points_array)
            if W is None:
                error_handling()
            print_matrix(W)
        case _:
            error_handling()