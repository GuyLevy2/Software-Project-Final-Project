import mykmeanssp as mksp
import numpy as np
import pandas as pd
import sys

def main():
    try:
        processedArgs = validateAndProcessInput(sys.argv)
        
        if len(processedArgs) == 0:
            print("Invalid Input!")
            return
        
        K = processedArgs[0]
        dimension = processedArgs[1]
        N = processedArgs[2]
        inputFileName = processedArgs[3]
        fileContent = processedArgs[4]
        goal = processedArgs[5]

        if goal == "spk":
            if K == 0:  # computes K by the eigengap heuristic
                K = mksp.eigengapHeuristic_fit(N, dimension, fileContent)
            
            # K-means++ algorithem
            returnValue_K_meansPP = k_meansPP(fileContent, N, K, dimension)
            if returnValue_K_meansPP == None:
                print("An Error Has Occurred")
                return
            
            centroids_keys = returnValue_K_meansPP[0]
            centroids_list = returnValue_K_meansPP[1]
            # OUTPUT - Print
            print(','.join(str(item) for item in centroids_keys))
            
            for i in range(K):
                print(','.join(str("%.4f"%item) for item in centroids_list[i]))

        elif goal == "wam":
            wam = mksp.wam_fit(N, dimension, fileContent)
            if wam == None:
                print("An Error Has Occurred")
                return

            for i in range(N):
                 print(','.join(str("%.4f"%item) for item in wam[i]))

        elif goal == "ddg":
            ddg = mksp.ddg_fit(N, dimension, fileContent)
            if ddg == None:
                print("An Error Has Occurred")
                return

            for i in range(N):
                 print(','.join(str("%.4f"%item) for item in ddg[i]))
        
        elif goal == "lnorm":
            lnorm = mksp.lnorm_fit(N, dimension, fileContent)
            if lnorm == None: 
                print("An Error Has Occurred")
                return
            
            for i in range(N):
                 print(','.join(str("%.4f"%item) for item in lnorm[i]))

        elif goal == "jacobi":
            eigenValues, eigenVectors = mksp.jacobi_fit(N, fileContent)
            # TODO check if returned value is None
            # TODO check the printing disign (eigen vectors should be printed as columns)

            print(','.join(str("%.4f"%item) for item in eigenValues))
            for i in range(N):
                print(','.join(str("%.4f"%item) for item in eigenVectors[i]))
    
    except Exception:
        print("An Error Has Occurred")
        return 
    
    return 

def k_meansPP(vectorsList, N, K, dimension):
    MAX_ITER = 300
    EPSILON = 0 # ???????TODO???????
    
    sortedOriginalKyes = vectorsList[:,0]    # The first column of vectors list represanted the original keys ordered

    np.random.seed(0) # SEED

    centroids_indexes = []
    centroids_keys = []

    indexes = range(N)
    D = [0.0 for i in indexes]

    i = 0
    random_centroid_ind = np.random.choice(N)
    centroids_indexes.append(random_centroid_ind)
    centroids_keys.append((int)(sortedOriginalKyes[random_centroid_ind]))

    while(i < K - 1):
        for ell in indexes:
            min_dist = -1

            for j in range(i+1):
                dist = np.square(np.linalg.norm(vectorsList[ell][1:] - vectorsList[centroids_indexes[j]][1:]))
                if min_dist == -1 or dist < min_dist:
                    min_dist = dist
            D[ell] = min_dist

        P = D / (np.sum(D)) # Probebilities array

        i+=1
        
        random_centroid_ind = np.random.choice(N, p=P)
        centroids_indexes.append(random_centroid_ind)
        centroids_keys.append((int)(sortedOriginalKyes[random_centroid_ind]))

    # Reordering vectorsList s.t. the first k vectors are the k cetroids have chosen in kmeans++ algorithem
    first_k_vectors = vectorsList[centroids_indexes]
    other_vectors = vectorsList[[i for i in indexes if i not in centroids_indexes]]
    vectorsList = np.concatenate((first_k_vectors, other_vectors), axis=0)

    # Cut the key value of vectorsList and convert numpu array into python list
    vectorsList = np.array(vectorsList)
    vectorsList = vectorsList[:,1:]
    vectorsList = vectorsList.tolist()

    # run kmeans algorithen throw C-API
    centroids_list = mksp.kmeans_fit(K, dimension, N, MAX_ITER, EPSILON, vectorsList)
    if centroids_list == None:
        return None # Rise error     
    
    return [centroids_keys, centroids_list]

def validateAndProcessInput(argsList):
    falseList = [] # will be returned in case of invalidation
    legalFileSuff = [".txt",".csv"]
    goalsList = ["spk", "wam", "ddg", "lnorm", "jacobi"]

    if len(argsList) != 4:
        return falseList

    k = str(argsList[1])
    goal = str(argsList[2])
    inputFileName = str(argsList[3])
    
    if goal not in goalsList:
        return falseList

    if len(inputFileName) < 5 or inputFileName[-4:] not in legalFileSuff:
        return falseList
    
    try:
        k = int(k) # is integer

        file = open(inputFileName,'rb')
        fileContent = pd.read_csv(file, header=None)
        fileContent = fileContent.to_numpy()
        
        #### Sotring the vectors by the column key ###
        #fileContent = fileContent[fileContent[:,0].argsort(),:] TODO not need
        dimension = fileContent.shape[1] - 1
        line_count = fileContent.shape[0]

        if k >= line_count or k < 0: # if K >= N
            return falseList


    except FileNotFoundError:
        return falseList

    except ValueError:
        return falseList
    
    except IndexError:
        return falseList

    except NameError:
        return falseList
    
    except Exception:
        return falseList
    
    return [k, dimension, line_count, inputFileName, fileContent, goal]

main()