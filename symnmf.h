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

void freeMatrix(float** matrix,int n);
float** matrixMultiplication(float** a, float** b, int n);
double distance(vector v1, vector v2);
void printMatrix(float** matrix, int n);
void freeMemory(all_vecs *all_vectors, int K, int N);
void printVector(vector *vec);
void errorHandling();
all_vecs getInput();
float **similarityMatrix(all_vecs points);
float **diagonalDegreeMatrix(all_vecs points);
float **normalizedSimilarityMatrix(all_vecs points);