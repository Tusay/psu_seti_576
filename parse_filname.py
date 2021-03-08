# This file contains a list of functions that can be used to parse and sort filterbank filenames

import numpy as np

#v This function takes in a filename, splits it into an array around the "." character and 
#v returns the two digit number in the first slot as a string:
def get_file_id(filename):
    return str(filename.split(".")[0][-2:])


#v This function takes in a filename, splits it into an array around the "_" character and 
#v returns the two digit number in the first slot as a string:
def get_blc(filename):
    return str(filename.split("_")[0][-2:])


def get_file_timestamp(filename):
    return str(''.join(filename.split("_")[-4:-2]))


def get_target_name(filename):
    return str(''.join(filename.split("_")[-2]))


#v This function sorts a list of files by the two digit file_id obtained in the get_file_id function above:
def sort_by_file_id(file_list):
    ids = []                                    #< Create an empty array
    for file in file_list:                      #< Start a for loop to iterate through the file_list
        ids.append(get_file_id(file))           #< Add the file_id of each file to the array of ids[]
    idx = np.argsort(ids)                       #< Sort the ids array of file_ids
    file_list = np.array(file_list)[idx]        #< Match the file_list array with the sorted idx array
    file_list = np.ndarray.tolist(file_list)    #< Update the file_list based on the sorted file_id array idx
    return file_list                            #< Return the sorted file_list


def sort_by_timestamp(file_list):
    times = []
    for file in file_list:
        times.append(get_file_timestamp(file))
    timex = np.argsort(times)
    file_list = np.array(file_list)[timex]
    file_list = np.ndarray.tolist(file_list)
    return file_list


#v This function sorts a list of files by the two digit blc number obtained in the get_blc function above:
def sort_by_blc(file_list):
    ids = []                                    # Same as the previous function
    for file in file_list:                      
        ids.append(get_blc(file))               
    idx = np.argsort(ids)                       
    file_list = np.array(file_list)[idx]    
    file_list = np.ndarray.tolist(file_list)    
    return file_list                            


def find_first_ON(dat_list, target_name):
    tar_list=[]
    for file in dat_list:
        tar_list.append(get_target_name(file))
    first_ON = tar_list.index(target_name)
    return first_ON
