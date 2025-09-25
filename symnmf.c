#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.h"

int initializeVector(int dim, vector *vec)
{
    vec->dimension = dim;
    vec->coordinates = (double *)calloc(dim, sizeof(double));
    if (vec->coordinates == NULL)
        return FAIL;
    return SUCCESS;
}

void freeVector(vector vec)
{
    if (vec.coordinates != NULL)
        free(vec.coordinates);
}

int initializeMatrix(int rows, int cols, matrix *newMatrix)
{
    int i;
    newMatrix->numOfRows = rows;
    newMatrix->numOfCols = cols;
    newMatrix->matrixEntries = (double **)calloc(newMatrix->numOfRows, sizeof(double *));
    if (newMatrix->matrixEntries == NULL)
        return FAIL;
    for (i = 0; i < newMatrix->numOfRows; i++)
    {
        newMatrix->matrixEntries[i] = (double *)calloc(newMatrix->numOfCols, sizeof(double));
        if (newMatrix->matrixEntries[i] == NULL)
        {
            freeMatrix(*newMatrix);
            return FAIL;
        }
    }
    return SUCCESS;
}

void freeMatrix(matrix mat)
{
    int i;
    if (mat.matrixEntries != NULL)
    {
        for (i = 0; i < mat.numOfRows; i++)
        {
            if (mat.matrixEntries[i] != NULL)
            {
                free(mat.matrixEntries[i]);
            }
        }
        free(mat.matrixEntries);
    }
}

int initializeDataPoints(int numOfVectors, int dim, dataPoints *points)
{
    int i;
    vector new_vec;
    points->numOfPoints = numOfVectors;
    points->points = (vector *)malloc(sizeof(vector) * numOfVectors);
    if (points->points == NULL)
        return FAIL;
    for (i = 0; i < points->numOfPoints; i++)
    {
        if (initializeVector(dim, &new_vec) == FAIL)
        {
            freeDataPoints(*points);
            return FAIL;
        }
        points->points[i] = new_vec;
    }
    return SUCCESS;
}

void freeDataPoints(dataPoints points)
{
    if (points.points != NULL)
    {
        int i;
        for (i = 0; i < points.numOfPoints; i++)
        {
            if (&points.points[i] != NULL)
                freeVector(points.points[i]);
        }
        free(points.points);
    }
}

/*matrixMultiplication calculates a*b, and saves the result in resMatrix*/
int matrixMultiplication(matrix a, matrix b, matrix *resMatrix)
{
    int i, j, k;
    double sum;
    if (a.numOfCols != b.numOfRows)
        return FAIL;
    if (initializeMatrix(a.numOfRows, b.numOfCols, resMatrix) == FAIL)
        return FAIL;
    for (i = 0; i < a.numOfRows; i++)
    {
        for (j = 0; j < b.numOfCols; j++)
        {
            sum = 0;
            for (k = 0; k < b.numOfRows; k++)
            {
                sum += a.matrixEntries[i][k] * b.matrixEntries[k][j];
            }
            resMatrix->matrixEntries[i][j] = sum;
        }
    }
    return SUCCESS;
}

/*distance calculates euclidian distance between two vectors*/
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

void printMatrix(matrix mat)
{
    int i, j;
    for (i = 0; i < mat.numOfRows; i++)
    {
        for (j = 0; j < mat.numOfCols; j++)
        {
            printf("%.4f", mat.matrixEntries[i][j]);
            if (j != mat.numOfCols - 1)
                printf(",");
            else
                printf("\n");
        }
    }
}

void printVector(vector vec)
{
    int i;
    for (i = 0; i < vec.dimension; i++)
    {
        if (i == vec.dimension - 1)
        {
            printf("%.4f\n", (vec.coordinates)[i]);
        }
        else
        {
            printf("%.4f,", (vec.coordinates)[i]);
        }
    }
}

void printDataPoints(dataPoints points)
{
    int i;
    for (i = 0; i < points.numOfPoints; i++)
    {
        printVector(points.points[i]);
    }
}

/*transpose computes transposed matrix mat, and saves the result in transposedMatrix*/
int transpose(matrix mat, matrix *transposedMatrix)
{
    int i, j;
    if (initializeMatrix(mat.numOfCols, mat.numOfRows, transposedMatrix) == FAIL)
        return FAIL;
    for (i = 0; i < mat.numOfCols; i++)
    {
        for (j = 0; j < mat.numOfRows; j++)
        {
            transposedMatrix->matrixEntries[i][j] = mat.matrixEntries[j][i];
        }
    }
    return SUCCESS;
}

/*trace calculates trace of matrix mat*/
double trace(matrix mat)
{
    double sum = 0;
    int i;
    for (i = 0; i < mat.numOfRows; i++)
    {
        sum += mat.matrixEntries[i][i];
    }
    return sum;
}

/*substractMatrices calculates A-B, and saves the result in result matrix*/
int substractMatrices(matrix A, matrix B, matrix *result)
{
    int i, j;
    if (A.numOfRows != B.numOfRows || A.numOfCols != B.numOfCols)
        return FAIL;
    if (initializeMatrix(A.numOfRows, A.numOfCols, result) == FAIL)
        return FAIL;
    for (i = 0; i < A.numOfRows; i++)
    {
        for (j = 0; j < A.numOfCols; j++)
        {
            result->matrixEntries[i][j] = A.matrixEntries[i][j] - B.matrixEntries[i][j];
        }
    }
    return SUCCESS;
}

/*updateH updates H as described in 1.4.2*/
int updateH(matrix H, matrix W, matrix *updatedH)
{
    matrix WmulH, transposedH, HMulTransposedH, HMulTransposedHMulH;
    int i, j;
    int result = SUCCESS;
    double denom;
    if (initializeMatrix(H.numOfRows, H.numOfCols, updatedH) == SUCCESS &&
        matrixMultiplication(W, H, &WmulH) == SUCCESS &&
        transpose(H, &transposedH) == SUCCESS &&
        matrixMultiplication(H, transposedH, &HMulTransposedH) == SUCCESS &&
        matrixMultiplication(HMulTransposedH, H, &HMulTransposedHMulH) == SUCCESS)
    {
        for (i = 0; i < H.numOfRows; i++)
        {
            for (j = 0; j < H.numOfCols; j++)
            {
                denom = HMulTransposedHMulH.matrixEntries[i][j];
                if (denom == 0)
                    denom += 0.000001;
                updatedH->matrixEntries[i][j] = H.matrixEntries[i][j] * (1 - BETA + BETA * (WmulH.matrixEntries[i][j] / denom));
            }
        }
    }
    else
        result = FAIL;
    freeMatrix(WmulH);
    freeMatrix(transposedH);
    freeMatrix(HMulTransposedH);
    freeMatrix(HMulTransposedHMulH);
    return result;
}

/*iterateAlgorithm updates H until max iteration number or convergence is reached*/
int iterateAlgorithm(matrix *H, matrix W)
{
    int i, j, l;
    matrix updatedH, updatedHMinusH, prevH;
    double frobeniusNorm;
    for (i = 0; i < MAXITER; i++)
    {
        if (updateH(*H, W, &updatedH) == FAIL)
            return FAIL;
        if (substractMatrices(updatedH, *H, &updatedHMinusH) == FAIL)
        {
            freeMatrix(updatedH);
            return FAIL;
        }
        frobeniusNorm = 0;
        for (j = 0; j < H->numOfRows; j++)
        {
            for (l = 0; l < H->numOfCols; l++)
            {
                frobeniusNorm += pow(fabs(updatedHMinusH.matrixEntries[j][l]), 2);
            }
        }
        freeMatrix(updatedHMinusH);
        prevH = *H;
        *H = updatedH;
        freeMatrix(prevH);
        if (frobeniusNorm < EPS)
            break;
    }
    return SUCCESS;
}

void errorHandling()
{
    printf("An Error Has Occurred\n");
    exit(1);
}

/*calculateNumOfPoints calculates the number of points in the given file*/
int calculateNumOfPoints(char *filename)
{
    int numOfPoints = 0;
    double n;
    char c;
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        return FAIL;
    }
    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
            numOfPoints++;
    }
    fclose(fp);
    return numOfPoints;
}

/*calculateDimension calculates the dimension of points in the given file*/
int calculateDimension(char *filename)
{
    int dimension = 0;
    FILE *fp = fopen(filename, "r");
    double n;
    char c;
    if (!fp)
    {
        return FAIL;
    }
    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        dimension++;
        if (c == '\n')
            break;
    }
    fclose(fp);
    return dimension;
}

int readPointsFromFile(char *filename, dataPoints *points)
{
    double coordinate;
    char c;
    int numOfPoints = calculateNumOfPoints(filename);
    int dimension = calculateDimension(filename);
    int i = 0, j = 0;
    FILE *fp = fopen(filename, "r");
    if (numOfPoints == FAIL)
        return FAIL;
    if (dimension == FAIL)
        return FAIL;
    if (initializeDataPoints(numOfPoints, dimension, points) == FAIL)
        return FAIL;
    if (!fp)
        return FAIL;
    while (fscanf(fp, "%lf%c", &coordinate, &c) == 2)
    {
        points->points[i].coordinates[j] = coordinate;
        j++;
        if (c == '\n')
        {
            i++;
            j = 0;
        }
    }
    fclose(fp);
    return SUCCESS;
}

/*similarityMatrix computes the similarity matrix and saves in in outputMatrix*/
int similarityMatrix(dataPoints points, matrix *outputMatrix)
{
    int i, j, n = points.numOfPoints;
    if (initializeMatrix(n, n, outputMatrix) == FAIL)
        return FAIL;
    for (i = 0; i < outputMatrix->numOfRows; i++)
    {
        for (j = 0; j < outputMatrix->numOfCols; j++)
        {
            if (j != i)
            {
                double exponent = -pow(distance(points.points[i], points.points[j]), 2) / 2;
                outputMatrix->matrixEntries[i][j] = exp(exponent);
            }
            else
                outputMatrix->matrixEntries[i][j] = 0;
        }
    }
    return SUCCESS;
}

/*diagonalDegreeMatrix calculates diagonal degree matrix of the given points and saves it in outputMatrix*/
int diagonalDegreeMatrix(dataPoints points, matrix *outputMatrix)
{
    int n = points.numOfPoints, i, j;
    matrix similarityMat;
    double rowSum;
    if (initializeMatrix(n, n, outputMatrix) == FAIL)
        return FAIL;
    if (similarityMatrix(points, &similarityMat) == FAIL)
    {
        freeMatrix(*outputMatrix);
        return FAIL;
    }
    for (i = 0; i < similarityMat.numOfRows; i++)
    {
        rowSum = 0;
        for (j = 0; j < similarityMat.numOfCols; j++)
        {
            rowSum += similarityMat.matrixEntries[i][j];
        }
        outputMatrix->matrixEntries[i][i] = rowSum;
    }
    freeMatrix(similarityMat);
    return SUCCESS;
}

/*normalizedSimilarityMatrix computes normalized similarity matrix of given points and saves the output in resMatrix*/
int normalizedSimilarityMatrix(dataPoints points, matrix *resMatrix)
{
    matrix similarityMat, degreeMat, mulMat;
    int i, result = SUCCESS;
    if (similarityMatrix(points, &similarityMat) == SUCCESS && diagonalDegreeMatrix(points, &degreeMat) == SUCCESS)
    {
        for (i = 0; i < degreeMat.numOfRows; i++)
        {
            degreeMat.matrixEntries[i][i] = 1 / (sqrt(degreeMat.matrixEntries[i][i]));
        }
        if (matrixMultiplication(degreeMat, similarityMat, &mulMat) == FAIL || matrixMultiplication(mulMat, degreeMat, resMatrix) == FAIL)
            result = FAIL;
    }
    else
        result = FAIL;
    freeMatrix(mulMat);
    freeMatrix(similarityMat);
    freeMatrix(degreeMat);
    return result;
}

int main(int argc, char **argv)
{
    char *goal;
    char *path;
    dataPoints points;
    matrix outputMatrix;
    if (argc != 3)
        return (1);
    goal = argv[1];
    path = argv[2];
    if (readPointsFromFile(path, &points) == FAIL)
        errorHandling();

    if (strcmp(goal, "sym") == 0)
    {
        if (similarityMatrix(points, &outputMatrix) == FAIL)
            errorHandling();
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        if (diagonalDegreeMatrix(points, &outputMatrix) == FAIL)
            errorHandling();
    }
    else if (strcmp(goal, "norm") == 0)
    {
        if (normalizedSimilarityMatrix(points, &outputMatrix) == FAIL)
            errorHandling();
    }
    else
        errorHandling();
    printMatrix(outputMatrix);
    freeMatrix(outputMatrix);
    freeDataPoints(points);
    return (0);
}