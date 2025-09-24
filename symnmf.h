typedef struct
{
    double *coordinates;
    int dimension;
} vector;

typedef struct
{
    /* TODO: change field name to 'vectors' ?*/
    vector *all_vectors;
    int num_vectors;
} dataPoints;
/* TODO: maybe change name to 'vectorList'? */

typedef struct
{
    double **matrixEntries;
    int numOfRows;
    int numOfCols;
} matrix;

int initializeVector(int dim, vector * vec);
void freeVector(vector vec);
int initializeDataPoints(int numOfVectors, int dim, dataPoints* points);
void freeDataPoints(dataPoints all_vectors);
int initializeMatrix(int rows, int cols, matrix *newMatrix);
void freeMatrix(matrix mat);
int matrixMultiplication(matrix a, matrix b, matrix *resMatrix);
double distance(vector v1, vector v2);
void printMatrix(matrix mat);
void printVector(vector vec);
void errorHandling();
int getInput(char *filename, dataPoints *points);
int similarityMatrix(dataPoints points, matrix *outputMatrix);
int diagonalDegreeMatrix(dataPoints points, matrix *outputMatrix);
int normalizedSimilarityMatrix(dataPoints points, matrix *resMatrix);
int transpose(matrix mat, matrix *transposedMatrix);
double trace(matrix mat);
int substractMatrices(matrix A, matrix B, matrix *result);
int updateH(matrix H, matrix W, matrix *updatedH);
int iterateAlgorithm(matrix *H, matrix W);
int calculateNumOfPoints(char *filename);
int calculateDimension(char *filename);
void printDataPoints(dataPoints points);

/* TODO: remove these when you are handing in the assignment? */
void testMatrixMultiplication();
void testDistance();
void testSimilarityMatrix();
void testDiagonalDegreeMatrix();
void testNormalizedSimilarityMatrix();
void testTranspose();
void testTrace();
void testSubstractMatrices();
void testUpdateH();