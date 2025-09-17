#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

void freeMatrix(double **matrix, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

double **matrixMultiplication(double **a, double **b, int rows1, int rows2, int cols2)
{
    double **resMatrix = (double **)calloc(rows1, sizeof(double *));
    int i;
    int j;
    for (i = 0; i < rows1; i++)
    {
        resMatrix[i] = (double *)calloc(cols2, sizeof(double));
        for (j = 0; j < cols2; j++)
        {
            int sum = 0;
            int k;
            for (k = 0; k < rows2; k++)
            {
                sum += a[i][k] * b[k][j];
            }
            resMatrix[i][j] = sum;
        }
    }
    return resMatrix;
}

double distance(vector v1, vector v2)
{
    double sum = 0;
    int i;
    for (i = 0; i < v1.dimension; i++)
    {
        sum += pow(v1.coordinates[i] - v2.coordinates[i], 2);
    }
    return sqrt(sum);
}

void printMatrix(double **matrix, int n, int k)
{
    int i;
    int j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            printf("%.4f", matrix[i][j]);
            if (j != k - 1)
                printf(",");
            else printf("\n");
        }
    }
}

void freeMemory(all_vecs *all_vectors, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        free(all_vectors->all_vectors[i].coordinates);
    }
    free(all_vectors->all_vectors);
}

void printVector(vector *vec)
{
    int i;
    for (i = 0; i < vec->dimension; i++)
    {
        if (i == vec->dimension - 1)
        {
            printf("%.4f", (vec->coordinates)[i]);
        }
        else
            printf("%.4f,", (vec->coordinates)[i]);
    }
    printf("\n");
}

void errorHandling()
{
    printf("An Error Has Occured\n");
}

all_vecs getInput()
{
    double n;
    char c;
    int i = 0, j = 0;
    all_vecs all_vectors;
    vector curr_vector;
    curr_vector.dimension = 0;
    curr_vector.coordinates = (double *)malloc(sizeof(double));
    if (curr_vector.coordinates == NULL)
        errorHandling();
    all_vectors.all_vectors = (vector *)malloc(sizeof(vector));
    if (all_vectors.all_vectors == NULL)
        errorHandling();
    while (scanf("%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
        {
            vector new_vector;
            curr_vector.coordinates[j] = n;
            j++;
            if (i == 0)
                curr_vector.dimension++;
            all_vectors.all_vectors[i] = curr_vector;
            i++;
            all_vectors.all_vectors = (vector *)realloc(all_vectors.all_vectors, sizeof(vector) * (i + 1));
            if (all_vectors.all_vectors == NULL)
                errorHandling();
            new_vector.dimension = j;
            new_vector.coordinates = (double *)malloc(sizeof(double) * new_vector.dimension);
            if (new_vector.coordinates == NULL)
                errorHandling();
            curr_vector = new_vector;
            j = 0;
            continue;
        }
        curr_vector.coordinates[j] = n;
        j++;
        if (i == 0)
        {
            curr_vector.dimension++;
            curr_vector.coordinates = (double *)realloc(curr_vector.coordinates, sizeof(double) * (j + 1));
            if (curr_vector.coordinates == NULL)
                errorHandling();
        }
    }
    free(curr_vector.coordinates);
    all_vectors.num_vectors = i;
    return all_vectors;
}

double **similarityMatrix(all_vecs points)
{
    int n = points.num_vectors;
    double **outputMatrix = (double **)calloc(n, sizeof(double *));
    int i;
    int j;
    for (i = 0; i < n; i++)
    {
        outputMatrix[i] = (double *)calloc(n, sizeof(double));
        for (j = 0; j < n; j++)
        {
            double exponent = -pow(distance(points.all_vectors[i], points.all_vectors[j]), 2) / 2;
            outputMatrix[i][j] = exp(exponent);
        }
    }
    return outputMatrix;
}
double **diagonalDegreeMatrix(all_vecs points)
{
    int n = points.num_vectors;
    double **outputMatrix = (double **)calloc(n, sizeof(double *));
    double **similarityMat = similarityMatrix(points);
    int i;
    int j;
    int rowSum;
    for (i = 0; i < n; i++)
    {
        rowSum = 0;
        outputMatrix[i] = (double *)calloc(n, sizeof(double));
        for (j = 0; j < n; j++)
        {
            rowSum += similarityMat[i][j];
        }
        outputMatrix[i][i] = rowSum;
    }
    freeMatrix(similarityMat, n);
    return outputMatrix;
}
double **normalizedSimilarityMatrix(all_vecs points)
{
    int n = points.num_vectors;
    double **similarityMat = similarityMatrix(points);
    double **degreeMat = diagonalDegreeMatrix(points);
    int i;
    double **mulMat;
    double **resMatrix;
    for (i = 0; i < n; i++)
    {
        degreeMat[i][i] = 1 / sqrt(degreeMat[i][i]);
    }
    mulMat = matrixMultiplication(degreeMat, similarityMat, n, n, n);
    resMatrix = matrixMultiplication(mulMat, degreeMat, n, n, n);
    freeMatrix(mulMat, n);
    freeMatrix(similarityMat, n);
    freeMatrix(degreeMat, n);
    return resMatrix;
}

/*int main(int argc, char **argv)
{
    char *goal = argv[1];
    all_vecs points = getInput();
    int n = points.num_vectors;
    double **outputMatrix;
    if(argc>3)
        return(1);
    char *path = argv[2];
    if (!strcmp(goal, "sym"))
    {
        outputMatrix = similarityMatrix(points);
    }
    if (!strcmp(goal, "ddg"))
    {
        outputMatrix = diagonalDegreeMatrix(points);
    }
    if (!strcmp(goal, "norm"))
    {
        outputMatrix = normalizedSimilarityMatrix(points);
    }
    printMatrix(outputMatrix, n, n);
    freeMatrix(outputMatrix, n);
    freeMemory(&points,n);
    return(0);
}
*/

int main()
{
    /*testMatrixMultiplication();*/
    /*testDistance();*/
    testSimilarityMatrix();
    return (0);
}