int wam_func(double*** vectors_list, int N, int dim, double*** outputWamMat);

int ddg_func(int N, double*** wamMat, double*** outputDdgMat);

int lNorm_func(int N, double*** wamMat, double*** ddgMat, double*** outputLNormMat);

int jacobi_func(int N, double*** symMat, double*** eigenVectors, double** eigenValues);

int eigenGap(int N, double** eigenValues);

int eigenComp(const void* a, const void* b);

double vectorDist(double* v1, double* v2, int dim);

int minusSqrtD(int N, double*** D);

int identityMat(int N, double*** mat);

int buildRotMat(int N, double*** A, int i, int j, int* c_p, int* s_p, double*** P);

int find_ij_pivot(int N, double*** A, int* i_p, int* j_p);

int matMult(int N, double*** mat1, double*** mat2, double*** outputMat);

int matDup(int N, double*** origMat, double*** newMat);

int matTranspose(int N, double*** mat, double*** matT);

int computeA_tag(int N, int i, int j, double*** A, double c, double s, double*** A_tag);

int convergenceTest(int N, double epsilon, double*** mat1, double*** mat2);

double offCalc(int N, double*** mat);

double*** initMat(int N);

int freeMat(int N, double*** mat);

int validateAndProcessInput(int argc, char* argv[], int* k, int* dimension, int* line_count, int* maxIter, char** inputFile, char** outputFile, double*** vectorsList);

int isValidInteger(char* str);