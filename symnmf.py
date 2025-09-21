import math
import sys
import numpy as np
import pandas as pd
import _symnmf

np.random.seed(1234)

def init_H(W, K):
    W_as_array = np.array(W, dtype=np.float64)
    n = W_as_array.shape[0]
    m = W_as_array.mean()
    upper = 2*np.sqrt(m/K)
    H = np.random.uniform(0, upper, size=(n, K))
    return H

def printmat(A):
    for row in A:
        line = ",".join(f"{val:.4f}" for val in row)
        print(line)

def getInputVariables():
    K = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return K, goal, file_name

def getDataPoints(file_name):
    vectors = pd.read_csv(file_name, header=None)
    vectors.columns = [f"coordinate {i}" for i in range(vectors.shape[1])]
    vectors.index.name = "key"
    return vectors

def getFinalH():
    W = _symnmf.norm(points_array)
    H = init_H(W, K)
    H = _symnmf.symnmf(H, W)
    return H

if __name__ == "__main__":
    K, goal, file_name = getInputVariables()
    vectors = getDataPoints(file_name)
    points_array=vectors.to_numpy()
    if goal == "symnmf":
        H = getFinalH()
        printmat(H)
    elif goal == "sym":
        A = _symnmf.sym(points_array)
        printmat(A)
    elif goal == "ddg":
        D = _symnmf.ddg(points_array)
        printmat(D)
    else:
        W = _symnmf.norm(points_array)
        printmat(W)