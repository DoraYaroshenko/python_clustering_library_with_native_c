#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

/*checked and workes -> printMatrix works, matrixMultiplication works*/
void testMatrixMultiplication()
{
    double **A = (double **)malloc(sizeof(double *) * 3);
    double **B = (double **)malloc(sizeof(double *) * 3);
    double **res;
    int i;
    int j;
    for (i = 0; i < 3; i++)
    {
        A[i] = (double *)malloc(sizeof(double) * 3);
        for (j = 0; j < 3; j++)
        {
            A[i][j] = i * j + j;
        }
    }
    printMatrix(A, 3, 3);
    for (i = 0; i < 3; i++)
    {
        B[i] = (double *)malloc(sizeof(double) * 3);
        for (j = 0; j < 3; j++)
        {
            B[i][j] = 1;
        }
    }
    printMatrix(B, 3, 3);
    res = matrixMultiplication(A, B, 3, 3, 3);
    printMatrix(res, 3, 3);
    freeMatrix(A, 3);
    freeMatrix(B, 3);
    freeMatrix(res, 3);
}

/*checked and works -> distance works and print vector works*/
void testDistance()
{
    vector a;
    vector b;
    int i;
    a.coordinates = (double *)malloc(sizeof(double) * 3);
    b.coordinates = (double *)malloc(sizeof(double) * 3);
    a.dimension = 3;
    b.dimension = 3;
    for (i = 0; i < 3; i++)
    {
        a.coordinates[i] = i;
        b.coordinates[i] = 2*i;
    }
    printVector(&a);
    printVector(&b);
    printf("%f",distance(a,b));
    free(a.coordinates);
    free(b.coordinates);
}

void testSimilarityMatrix(){
    int i;
    int j;
    double** similarityMat;
    all_vecs points;
    points.all_vectors = (vector*)malloc(sizeof(vector)*3);
    points.num_vectors = 3;
    for(i=0;i<3;i++){
        points.all_vectors[i].dimension = 4;
        points.all_vectors[i].coordinates = (double*)malloc(sizeof(double)*4);
        for(j=0;j<4;j++){
            points.all_vectors[i].coordinates[j] = i+j;
        }
        printVector(&points.all_vectors[i]);
    }
    similarityMat = similarityMatrix(points);
    printMatrix(similarityMat,3,3);
    freeMemory(&points,3);
    freeMatrix(similarityMat,3);
}