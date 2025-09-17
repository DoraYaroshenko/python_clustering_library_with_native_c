import math
import sys
import numpy as np
import pandas as pd
import symnmf

np.random.seed(1234)

# def similarity_matrix(points_array):
#     squares_sum = np.sum(points_array**2, axis=1, keepdims=True) # לכל שורה במטריצה מעלה את האיבר בכל עמודה בריבוע וסוכם את העמודות. זה נהיה וקטור עמודה שמכיל את הסכומים הנ״ל בכל שורה
#     squared_sub_mat = squares_sum+squares_sum.T # זה מטריצה שבה xij = |x_i|^2+|x_j|^2
#     squared_sub_mat = squared_sub_mat-2*(points_array @ points_array.T) # מכל איבר איקס איי ג׳יי במטריצה הקודמת יצא לנו להחסיר שתיים כפול מ״פ בין איקס איי לאיקס ג׳יי 
#     # סהכ מתקיים מ״חוק הקוסינוסים״ ש*לא* למדנו עם סמיון: ||x_i - x_j||^2 = ||x_i|| + ||x_j|| - 2*<x_i, x_j> and it matches what we have done so far
#     A = np.exp(-squared_sub_mat/2) # מטריצה לפי ההגדרה, עד כדי זה שהאלכסון לא תואם להגדרה
#     np.fill_diagonal(A, 0) # איפוס האלכסון
#     return A

# def degree_matrix(A):
#     d = np.sum(A, axis=1) # סוכם את האיברים בכל שורה
#     D = np.diag(d) # מציב את האיברים של די, שזה הדרגות שחישבנו קודם, על האלכסון - ושאר המטריצה היא אפסים
#     return D

# def normalized_similarity_mat(D, A):
#     d = np.diag(D)
#     D_power = np.diag(1.0 / np.sqrt(d)) # מעלים את האלכסון בחזקת מינוס חצי שזה כמו אחד חלקי שורש ומציבים אותו באלכסון המטריצה החדשה שלנו
#     W = D_inv_sqrt @ A @ D_inv_sqrt # מחשבים את הביטוי שהוגדר
#     return W

def init_H(W, K):
    n = W.shape[0]
    m = W.mean()
    upper = 2*np.sqrt(m/K)
    H = np.random.uniform(0, upper, size=(n, K))
    return H

def printmat(A):
    for row in A:
        line = ",".join(f"{val:.4f}" for val in row)
        print(line)

def getInputVariables():
    #it has to be exactly 4 args - I hope we still need to check that? because apparently there is no need to check the arguments themselves in this assignment?
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)

    K = float(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return K, goal, file_name

def getDataPoints(file_name):
    vectors = pd.read_csv(file_name, header=None)
    vectors.columns = [f"coordinate {i}" for i in range(df.shape[1])]
    vectors.index.name = "key"
    return vectors

if _name_ == "_main_":
    K, goal, file_name = getInputVariables() #ok
    vectors = getDataPoints(file_name) #ok
    N = len(vectors) #ok
    points_array=vectors.to_numpy() # okay, returns actual matrix in numpy form..

    if goal == "symnmf":
        A = symnmf.sym(points_array) # maybe should use sym method from c file instead.
        D = symnmf.ddg(A) # maybe should use ddg method from c file instead
        W = symnmf.norm(D, A) # should have used norm from c file instead apparently
        H = init_H(W, K)
        H = symnmf.symnmf(H, W)
        printmat(H)
    elif goal == "sym":
        A = symnmf.sym(points_array)
        printmat(A)
    elif goal == "ddg":
        D = symnmf.ddg(points_array)
        printmat(D)
    else: # goal == "norm"
        W = symnmf.norm(points_array)
        printmat(W)