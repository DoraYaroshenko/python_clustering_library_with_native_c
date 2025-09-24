import math
import sys
import numpy as np
import pandas as pd
import _symnmf

np.random.seed(1234)

def errorHandling():
    print("An Error Has Occurred")
    sys.exit(1)

def init_H(W, K):
    W_as_array = np.array(W)
    n = W_as_array.shape[0]
    m = W_as_array.mean()
    upper = 2*np.sqrt(m/K)
    H = np.random.uniform(0, upper, size=(n, K))
    return H

# TODO: change name to 'printMatrix'
def printmat(A):
    for row in A:
        line = ",".join(f"{val:.4f}" for val in row)
        print(line)

def getInputVariables():
    # TODO: validate number of input variables?
    K = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return K, goal, file_name

def getDataPoints(file_name):
    vectors = pd.read_csv(file_name, header=None)
    vectors.columns = [f"coordinate {i}" for i in range(vectors.shape[1])]
    vectors.index.name = "key"
    return vectors

# what is final H?
# TODO: maybe change name to 'calculateDecompositionMatrix' or just 'calculateH'
def getFinalH(points_array, K):
    W = _symnmf.norm(points_array)
    if W is None:
        errorHandling()
    H = init_H(W, K)
    H = _symnmf.symnmf(H, W)
    if H is None:
        errorHandling()
    return H

if __name__ == "__main__":
    # general comment: in python the convention is snake case
    # get_input_variables
    # getInputVariables - is in java or smth, not python
    # you dont have to change it, this is just so you know
    K, goal, file_name = getInputVariables()
    vectors = getDataPoints(file_name)
    N = len(vectors)
    if K>=N or K<=1:
        errorHandling()
    
    points_array=vectors.to_numpy()
    match goal:
        case "sym":
            A = _symnmf.sym(points_array)
            if A is None:
                errorHandling()
            printmat(A)
        case "symnmf":
            H = getFinalH(points_array, K)
            printmat(H)
        case "ddg":
            D = _symnmf.ddg(points_array)
            if D is None:
                errorHandling()
            printmat(D)
        case "norm":
            W = _symnmf.norm(points_array)
            if W is None:
                errorHandling()
            printmat(W)
        case _:
            errorHandling()
    sys.exit(0) # remove - not needed