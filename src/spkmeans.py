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
        maxIter = processedArgs[3]
        epsilon = processedArgs[4]
        inputFileName1 = processedArgs[5]
        inputFileName2 = processedArgs[6]
        vectorsList = processedArgs[7]

        sortedOriginalKyes= vectorsList[:,0]    # The first column of vectors list represanted the original keys ordered
        
        # K-means++ algorithem
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
        centroids_list = mksp.fit(K, dimension, N, maxIter, epsilon, vectorsList) # check what happens if null is returned. TEST
        if centroids_list == None:
            print("An Error Has Occurred")
            return 

        
        # OUTPUT - Print
        print(','.join(str(item) for item in centroids_keys))
        
        for i in range(K):
            print(','.join(str("%.4f"%item) for item in centroids_list[i]))

    except Exception:
        print("An Error Has Occurred")
        return 
    
    return 

def validateAndProcessInput(argsList):
    falseList = [] # will be returned in case of invalidation

    if len(argsList) > 6 or len(argsList) < 5:
        return falseList

    if len(argsList) == 5:
        k = str(argsList[1])
        maxIter = str(300)
        epsilon = str(argsList[2])
        inputFileName1 = str(argsList[3])
        inputFileName2 = str(argsList[4])

    if len(argsList) == 6:
        k = str(argsList[1])
        maxIter = str(argsList[2])
        epsilon = str(argsList[3])
        inputFileName1 = str(argsList[4])
        inputFileName2 = str(argsList[5])

    
    if (len(inputFileName1) < 5 or inputFileName1[-4:] not in [".txt",".csv"]):
        return falseList
    
    if (len(inputFileName2) < 5 or inputFileName2[-4:] not in [".txt",".csv"]):
        return falseList
    
    try:
        k = int(k) # is integer
        maxIter = int(maxIter) # is integer
        epsilon = float(epsilon)

        if k <= 1 or maxIter < 0: # if k <= 1 and maxIter < 0
            return falseList
        
        if epsilon < 0:
            return falseList


        file1 = open(inputFileName1,'rb')
        vectorsList1 = pd.read_csv(file1, header=None)
        file2 = open(inputFileName2,'rb')
        vectorsList2 = pd.read_csv(file2, header=None)

        vectors_merge = pd.merge(vectorsList1, vectorsList2, how='inner', on=0)        
        vectorsList = vectors_merge.to_numpy()
        
        #### Sotring the vectors by the column key ###
        vectorsList = vectorsList[vectorsList[:,0].argsort(),:]
        dimension = vectorsList.shape[1] - 1
        line_count = vectorsList.shape[0]

        if k >= line_count: # if K >= N
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
    
    return [k, dimension, line_count ,maxIter, epsilon, inputFileName1, inputFileName2, vectorsList]

main()
