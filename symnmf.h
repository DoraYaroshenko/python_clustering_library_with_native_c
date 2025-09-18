typedef struct
{
    double *coordinates;
    int dimension;
} vector;

typedef struct
{
    vector *all_vectors;
    int num_vectors;
} all_vecs;

void freeMatrix(double** matrix,int n);
double **matrixMultiplication(double **a, double **b, int rows1, int rows2, int cols2);
double distance(vector v1, vector v2);
void printMatrix(double** matrix, int n, int k);
void freeMemory(all_vecs *all_vectors, int N);
void printVector(vector *vec);
void errorHandling();
all_vecs getInput();
double **similarityMatrix(all_vecs points);
double **diagonalDegreeMatrix(all_vecs points);
double **normalizedSimilarityMatrix(all_vecs points);
double **transpose(double **matrix, int rows, int cols);
double trace(double **matrix, int n);
double **substractMatrices(double **A, double **B, int rows, int cols);
double **updateH(double **H, double **W, int n, int k);
double **iterateAlgorithm(double **H, double **W, int n, int k);

void testMatrixMultiplication();
void testDistance();
void testSimilarityMatrix();
void testDiagonalDegreeMatrix();
void testNormalizedSimilarityMatrix();
void testTranspose();
void testTrace();
void testSubstractMatrices();
void testUpdateH();