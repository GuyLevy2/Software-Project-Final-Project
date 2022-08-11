#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
    int k = 0, dimension, line_count;
    int maxIter, iteration_number;
    char* inputFile, *outputFile;
    double** vectorsList, **cluster, **centroidsList;;
    double EPSILON = 0.001;
    int inputError, writeSuccess, freeSuccess;
    double*** clusterList;
    int* index_to_insert;
    int i, j;
    double* vec, *centroid_i;
    double deltaMiu;

    inputError = validateAndProcessInput(argc, argv, &k, &dimension, &line_count, &maxIter, &inputFile, &outputFile, &vectorsList);
    
    if (!inputError){
        printf("Invalid Input!");
        return 1;
    }
}

/* 
 * Function: validateAndProcessInput
 * ---------------------------------
 * processes the CMD input and file
 *  
 * argc: argc
 * argv: argv
 * k: a pointer to k
 * dimension: a pointer to the dimension
 * line_count: a pointer to N
 * maxIter: a pointer to the max iterations number
 * inputFile: a pointer to the input file path
 * outputfile: a pointer to the output file path
 * vectors_list: a pointer to the list of all N vectors
 * 
 * returns: 0 if the input is invalid and 1 if the input is valid
 */
int validateAndProcessInput(int argc, char* argv[], int* k, int* dimension, int* line_count, int* maxIter, char** inputFile, char** outputFile, double*** vectorsList){
    int isTxtFile; /* new */
    int isCSVFile; /* new */
    char *k_str;
    char *maxIter_str;
    int outputStrLen;
    char c;
    int after_firstline;
    int curr_dimension;
    double *vec;
    int vectors_index;
    int vectorList_index;
    int i, j;
    FILE *fp1, *fp2;
    int ch; /* Guy added because of fgetc func giving error EOF??? */

    if (argc > 5 || argc < 4){
        return 0;
    }

    if (argc == 4){
        k_str = argv[1];
        maxIter_str = "200";
        *inputFile = argv[2];
        *outputFile = argv[3];
    }

    if (argc == 5){
        k_str = argv[1];
        maxIter_str = argv[2];
        *inputFile = argv[3];
        *outputFile = argv[4];
    }

    if(!isValidInteger(k_str)){ /* checking that k_str is a valid number - may cause problems and might not be needed? */
        return 0;
    }
    
    *k = atoi(k_str);
    
    if(!isValidInteger(maxIter_str)){ /* checking that str is a valid number - is it needed? */
        return 0;
    }
    
    *maxIter = atoi(maxIter_str);

    if (*k <= 1 || *maxIter < 0){
        return 0;
    }

    /* checking that the output file is in the right format */
    outputStrLen = 0;
    while ((*outputFile)[outputStrLen] != '\0'){
        outputStrLen++;
    }

    if (outputStrLen < 5){
        return 0;
    }

    isTxtFile = 0;
    isCSVFile = 0;

    if ((*outputFile)[outputStrLen - 1] == 't' 
        || (*outputFile)[outputStrLen - 2] == 'x' 
        || (*outputFile)[outputStrLen - 3] == 't' 
        || (*outputFile)[outputStrLen - 4] == '.'
        ){
            isTxtFile = 1;
        }
    
    if ((*outputFile)[outputStrLen - 1] == 'c' 
        || (*outputFile)[outputStrLen - 2] == 's' 
        || (*outputFile)[outputStrLen - 3] == 'v' 
        || (*outputFile)[outputStrLen - 4] == '.'
        ){
            isCSVFile = 1;
        return 0;
    }

    if (isTxtFile == 0 && isCSVFile == 0){
        return 0;
    }
    
    fp1 = fopen(*inputFile, "r"); /* First open */
    if (fp1 == NULL){ /* if can't open the file. */
        return 0;
    }
    
    if (isTxtFile == 1){
        *dimension = 1;
        *line_count = 0;
        
        after_firstline = 0;
        curr_dimension = 1;
    
        for (ch = getc(fp1); ch != EOF; ch = getc(fp1)){ /* First run - dimension and number of vectors. */
            c = (char) ch; /* Guy added because of fgetc func giving error EOF - mem leak??? */
            if (c == ','){
                if (after_firstline == 0){
                    *dimension += 1;
                }
                else{
                    curr_dimension += 1;
                }
            }

            if (c == '\n'){
                *line_count = *line_count + 1;
                after_firstline = 1;
                
                if(*dimension != curr_dimension && curr_dimension != 1){
                    fclose(fp1); /* Guy Memory leak */
                    return 0;
                }
                
                curr_dimension = 1;
            }

            /* if there is a char that is not a number or , or \n. */
            if ((c != '.') && (c != '-') && (c != ',') && (c != '\n') && (isdigit(c) == 0)){
                fclose(fp1); /* `Guy Memory leak problem */
                return 0;
            }
        }
        
        fclose(fp1);
        
        if (*k >= *line_count){ /* if K >= N */
                return 0;
        } 

        fp2 = fopen(*inputFile, "r"); /* Second open */
        if (fp2 == NULL){ /* if can't open the file. */
            return 0;
        }

        *vectorsList = malloc(*line_count * sizeof(double*));
        vec = malloc(*dimension * sizeof(double));
        vectors_index = 0;
        vectorList_index = 0;
        
        for(i = 0; i < *line_count; i++){
            for (j = 0; j < *dimension - 1; j++){
                fscanf(fp2, "%lf,", &(vec[vectors_index]));
                vectors_index++;
            }

            fscanf(fp2, "%lf\n", &(vec[vectors_index]));
            vectors_index++;

            (*vectorsList)[vectorList_index] = vec;
            vectorList_index++;
            
            vectors_index = 0;
            vec = malloc(*dimension * sizeof(double));
        }

        free(vec); /* Guy Mem Testing */

        fclose(fp2);

    }

    return 1;
}

/* 
 * Function: isValidInteger
 * ------------------------
 * checks if a string represents a valid int
 *  
 * str: the string to be checked
 * 
 * returns: 1 if the string is a valid int, and 0 otherwise
 */
int isValidInteger(char* str){
    int i;
    char c;

    if (str[0] == '\0'){
        return 0;
    }
    
    i = 0;
    c = str[0];
    while (c != '\0'){ /* checking that str is a valid number */
      if (c < 48 || c > 57){
         return 0; /* str is not an integer */
      }
      else{
         i++;
      }
    c = str[i];
   }

   return 1;
}