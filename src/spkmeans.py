import spkmeansmodule as mksp
import numpy as np
import pandas as pd
import sys

def main():
    try:
        processedArgs = validateAndProcessInput(sys.argv)
        if len(processedArgs) == 0: # processedArgs = [] in case of invalid input
            print("Invalid Input!")
            return
        
        K = processedArgs[0]
        dimension = processedArgs[1]
        N = processedArgs[2]
        inputFileName = processedArgs[3]
        fileContent = processedArgs[4]
        goal = processedArgs[5]

        if goal == "spk":
            # The matrix T is the output after stages (1)-(5) of the spk algorithem
            T = mksp.spk_fit(N, dimension, K, fileContent)
            if T == None:
                raise Exception
            
            K = np.asarray(T).shape[1] # The final K in the number of columns of the matrix T

            # Error handling
            if K == 1:
                raise Exception 
                       
            # stage (6): Treating each row of T as a point in RxK, cluster them into k clusters via the K-means algorithm
            returnValue_K_meansPP = k_meansPP(T, N, K, K)
            
            # Error handling
            if returnValue_K_meansPP == None:
                raise Exception

            centroids_keys = returnValue_K_meansPP[0]
            centroids_list = returnValue_K_meansPP[1]

            # Output - printing
            print(','.join(str(item) for item in centroids_keys))
            printFloatMatrix_format(K, centroids_list)

        elif goal == "wam":
            wam = mksp.wam_fit(N, dimension, fileContent)
            
            # Error handling
            if wam == None:
                raise Exception
            
            # Output - printing
            printFloatMatrix_format(N, wam)

        elif goal == "ddg":
            ddg = mksp.ddg_fit(N, dimension, fileContent)
            
            # Error handling
            if ddg == None:
                raise Exception

            # Output - printing
            printFloatMatrix_format(N, ddg)
        
        elif goal == "lnorm":
            lnorm = mksp.lnorm_fit(N, dimension, fileContent)
            
            # Error handling
            if lnorm == None: 
                raise Exception
            
            # Output - printing
            printFloatMatrix_format(N, lnorm)

        elif goal == "jacobi":
            eigenValues, eigenVectors = mksp.jacobi_fit(N, fileContent)

            # Error handling
            if eigenValues == None or eigenVectors == None:
                raise Exception

            # Output - printing
            print(','.join(str("%.4f"%item) for item in eigenValues))
            printFloatMatrix_format(N, eigenVectors)
    
    except:
        print("An Error Has Occurred")
        return 
    
    return 

def printFloatMatrix_format(numOfRows, matrix):
    """ 
    Function: printFloatMatrix_format:
    ----------------------------------
    prints a matrix with float etries with a `,` between 2 following etries
        and `\n` between 2 following rows

    Parameters:
    -----------
    numOfRows: int
        The number of rows in matrix
    matrix: list
        A 2-dimentional matrix with float etries
    
    Return: None
    """
    for i in range(numOfRows):
        print(','.join(str("%.4f"%item) for item in matrix[i]))
    
    return

def k_meansPP(vectorsList, N, K, dimension):
    """ 
    Function: k_meansPP:
    --------------------
    runs k-means++ algorithm:
    1. Computes initial k cetroids for k-means algorithm
    2. Runs k-means algorithm

    Parameters:
    -----------
    vectorsList: list
        A list of N vectors, each of them with float etries
    N: int
        The number of vectors in vectorsList
    K: int
        The number of initial k-centroids
    dimention: int
        The dimention of the vectors
    
    Return:
    -------
    list
        a list of 2 lists:
        1. A list of the initial k centroids indexes
        2. A list of the final k centroids after k-means running
    returns None in case of error durnig the function
    """

    MAX_ITER = 300
    EPSILON = 0.0 
    np.random.seed(0) # SEED
    centroids_indexes = []  # index from 0 to N-1
    indexes = range(N)
    D = [0.0 for i in indexes]
    
    i = 0
    random_centroid_ind = np.random.choice(N)
    centroids_indexes.append(random_centroid_ind)

    # K means ++ iterations
    while(i < K - 1):
        # computing distance array
        for ell in indexes:
            min_dist = -1
            for j in range(i + 1):
                v_ell = np.asarray(vectorsList[ell])
                v_cent_j = np.asarray(vectorsList[centroids_indexes[j]])
                dist = np.square(np.linalg.norm(v_ell - v_cent_j))
                if min_dist == -1 or dist < min_dist:
                    min_dist = dist
            D[ell] = min_dist
        # computing Probebilities array
        P = D / (np.sum(D))
        i+=1
        # choosing centroid
        random_centroid_ind = np.random.choice(N, p=P)
        centroids_indexes.append(random_centroid_ind)

    # Reordering vectorsList s.t. the first k vectors are the k cetroids have chosen in kmeans++ algorithem
    vectorsList = np.array(vectorsList)
    first_k_vectors = vectorsList[centroids_indexes]
    other_vectors = vectorsList[[i for i in indexes if i not in centroids_indexes]]
    vectorsList = np.concatenate((first_k_vectors, other_vectors), axis=0)
    vectorsList = vectorsList.tolist()
    
    # run kmeans algorithen throw C-API
    centroids_list = mksp.kmeans_fit(K, dimension, N, MAX_ITER, EPSILON, vectorsList)
    if centroids_list == None:
        return None # Raise error     

    return [centroids_indexes, centroids_list]

def validateAndProcessInput(argsList):
    """ 
    Function: validateAndProcessInput:
    ----------------------------------
    validates the initial arguments of the program

    Parameters:
    -----------
    argsList: list
        A list of the program's arguments

    Return:
    -------
    list
        a list of 2 lists:
        1. A list of the initial k centroids indexes
        2. A list of the final k centroids after k-means running
    returns [] in case of invalid input.
    The input arguments must be 4 elments:
    1. The name of the program
    2. An integer-values parameter k
    3. The goal of the program from GOALSLIST
    4. A path to .txt or .csv file contains float-valued vectors
    """

    falseList = [] # will be returned in case of invalidation
    LEGALFILES = [".txt",".csv"]
    GOALSLIST = ["spk", "wam", "ddg", "lnorm", "jacobi"]

    if len(argsList) != 4:
        return falseList

    k = str(argsList[1])
    goal = str(argsList[2])
    inputFileName = str(argsList[3])
    
    if goal not in GOALSLIST:
        return falseList

    if len(inputFileName) < 5 or inputFileName[-4:] not in LEGALFILES:
        return falseList
    
    try:
        k = int(k) # is integer

        file = open(inputFileName,'rb')
        fileContent = pd.read_csv(file, header=None)
        fileContent = fileContent.to_numpy()

        dimension = fileContent.shape[1]
        line_count = fileContent.shape[0]

        if k >= line_count or k < 0:
            return falseList
  
    except:
        return falseList
    
    return [k, dimension, line_count, inputFileName, fileContent.tolist(), goal]

main()