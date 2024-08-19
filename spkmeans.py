import myspkmeans as spkmeans
import numpy as np
import pandas as pd
import argparse
import sys
import math

max_iter = 300 #global

def printList(list):
    for centroid in list:
         print(','.join([format(centroid[j], ".4f") for j in range(len(centroid))]))

def printDiagonal(list):
     print(','.join([format(list[j][j], ".4f") for j in range(len(list))]))

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

def print_required(list1, list2):
    for i in range(len( list1)):
        if (i != (len( list1) - 1)):
            print( list1[i], end=",")
        else:
            print( list1[i])

    printList(list2)

def distance_calculator(v1, v2):
    assert (len(v1) == len(v2))
    distance = 0
    for c in range(len(v1)):
        coordinates_d = (v1[c] - v2[c]) ** 2
        distance += coordinates_d

    return distance

def probability_calculator(i, D_lst, D_sum):
    prob = float ( (D_lst[i] / D_sum) )
    return prob

def kmeans_pp(k, data_array, n, d):

    np.random.seed(0)
    first_centroid_index = np.random.choice(n)  #choose row number randomally
    initial_centroids_indices = [-1 for i in range(k+1)]  
    initial_centroids_indices[1] = first_centroid_index
    D_lst = [0 for i in range(n)]
    P_lst = [0 for i in range(n)]
    z = 1

    first_k_centroids = np.empty((0, d))
    first_k_centroids = np.append(first_k_centroids, [data_array[first_centroid_index]], axis=0)


    while ( z < k ):
            D_sum = 0
            for i in range(n):   #for each xi 
                min_distance = math.inf

                for j in range(1,z+1):   
                        v1 = data_array[i]
                        indice_of_chosen_vector_in_j_th_folder = initial_centroids_indices[j]
                        v2 = data_array[indice_of_chosen_vector_in_j_th_folder]
                        curr_distance = distance_calculator(v1,v2)

                        if(curr_distance<=min_distance):
                            min_distance=curr_distance

                D_lst[i] = min_distance

            for i in range(n):   #calculate D_sum
                D_sum += D_lst[i]
            for i in range(n):  #find P(Xi) for each Xi and update P_lst
                P_lst[i] = probability_calculator(i,D_lst,D_sum)

            z += 1
            z_centroid_index = np.random.choice(n, p=P_lst)
            initial_centroids_indices[z] = z_centroid_index
            first_k_centroids = np.append(first_k_centroids, [data_array[z_centroid_index]], axis=0)

    return first_k_centroids, initial_centroids_indices[1:len(initial_centroids_indices)]
    
def receiving_user_inputs():
    if len(sys.argv) != 4:
        print("Invalid Input!")
        return
    
    reader = argparse.ArgumentParser()
    reader.add_argument("k", type = int)
    reader.add_argument("goal", type = str)
    reader.add_argument("fileName", type = str)
    arguments = reader.parse_args()
    return arguments

def mainFunc():
    arguments = receiving_user_inputs()
    if (arguments.k < 0):
        print("Invalid Input")
        return
    if arguments.goal not in {"spk", "lnorm", "wam", "ddg", "jacobi"}:
        print("Invalid Input")
        return

    db = pd.read_csv(arguments.fileName, header=None)
    pyList = db.to_numpy().tolist()
    if arguments.goal == "wam":
        res = spkmeans.calcWAM(pyList, len(pyList), len(pyList[0]))
        printList(res)
    if arguments.goal == "ddg":
        if (len(pyList) != len(pyList[0])):
            print("Invalid Input")
            return
        res = spkmeans.calcDDG(pyList, len(pyList), len(pyList[0])) #need further testing
        printList(res)
    if arguments.goal == "lnorm":
        res = spkmeans.calcLnorm(pyList, len(pyList), len(pyList[0]))
        printList(res)
    if arguments.goal == "jacobi":
        if len(pyList) != len(pyList[0]) or not check_symmetric(db.to_numpy()):
            print("Invalid Input")
            return
        res = spkmeans.calcJacobi(pyList, len(pyList), len(pyList[0]))
        print("eigenValues:")
        printDiagonal(res[0])
        print("eigenVectors:")
        printList(res[1])
    if arguments.goal == "spk":
        if arguments.k == 0:
            inpMAt = spkmeans.calcTMatrix(pyList, len(pyList), len(pyList[0]), -1)
            numOfClusters = len(inpMAt[0])
        else:
            numOfClusters = arguments.k
            inpMAt = spkmeans.calcTMatrix(pyList, len(pyList), len(pyList[0]), numOfClusters)
        numOfVectors = len(inpMAt)
        first_k_centroids, first_k_centroids_indices =  kmeans_pp(numOfClusters, inpMAt, numOfVectors, len(inpMAt[0]))
        first_k_centroids = first_k_centroids.tolist()
        the_final_centroids = spkmeans.kmeansppCalc(first_k_centroids, inpMAt, max_iter, numOfVectors, len(inpMAt[0]), numOfClusters, 0.005) # here we do the actual calling

        if (the_final_centroids == None):
            print("An Error Has Occured")
        print_required(first_k_centroids_indices, the_final_centroids)        

mainFunc()