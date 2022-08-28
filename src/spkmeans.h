#ifndef SPKMEANS_H
#define SPKMEANS_H

int wam_func(double***, int, int, double***);
int ddg_func(int, double***, double***);
int lNorm_func(int, double***, double***, double***);
int jacobi_func(int, double***, double***, double**);
int kmeans_c(int, int, int, int, double, double**, double***);

double** initMat(int);
double** initMatMN(int, int);
void freeMat(int, double***);

int eigenGap(int, double**);
int sortEigenValuesAndEigenVectors(int, double**, double***);
int Fill_K_LargestEigenVectors(int, int, double***, double***);
int ReNormalizedRows(int, int, double***, double***);

#endif