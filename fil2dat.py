# This code is designed to turn filterbank .fil files into .dat files using turboSETI.
# It will output .h5, .dat, and .log files into an output '/processed' folder as well as a list of .fil files used in a fil_files.lst file.
#
# The input arguments for the code are currently:
#       1. '/input_directory/' where all the filterbank .fil files are located, with or without the trailing slash.
#       2. '--clobber' if present this will purge the output directory of all files. 
#          Otherwise, it will check if the files are already present in the output directory and skip over them if so.
#       3. '--gpu' utilization is technically possible with turboSETI, though untested by myself as of yet. Will test and update this later.

#v PACKAGES:
import numpy as np
import os
import sys
import glob
import time
import argparse
import csv
import logging
import pandas as pd
import turbo_seti
import blimpy as bp
from subprocess import call
from drift_rate import get_drift_rate
from parse_filname import get_file_id
from parse_filname import sort_by_timestamp
from turbo_seti.find_doppler.find_doppler import FindDoppler

#v The pickle package and function are for code testing purposes only:
import pickle
def my_little_pickle(thingy, name):
    print(f"\nMy little pickle: {name}\n{pickle.loads(pickle.dumps(thingy))}\n") 

#v FUNCTIONS:
#v This function parses the input arguments:
def parse_args():
    parser = argparse.ArgumentParser(description='Process BL filterbank data.')
    parser.add_argument('indir', metavar='/input_data_directory/', type=str, nargs=1,
                        help='directory containing the .fil files')
    parser.add_argument('-clobber', action='store_true',
                        help='overwrite files if they already exist')
    parser.add_argument('-GPU', action='store_true',
                        help='use GPU acceleration if possible')
    parser.add_argument('-max', dest='MaxDF', metavar='#Number', type=float, default=1,
                        help='set the max drift rate factor as number')
    parser.add_argument('-min', dest='MinDF', metavar='#Number', type=float, default=0,
                        help='set the minimum drift rate factor as number')
    parser.add_argument('-s1', dest='SNR1', metavar='#Number', type=float, default=10,
                        help='set the signal to noise ratio as a number to search for hits')
    parser.add_argument('-out', metavar='/output_folder/', type=str, default='processed/', nargs=1,
                        help='folder within indir to put output files into')
    args = parser.parse_args()
    # Check for trailing slash in the directory path and add it if absent
    odict = vars(args)
    if odict["indir"]:
        indir = odict["indir"][0]
        if indir[-1] != "/":
            indir += "/"
        odict["indir"] = indir
    if odict["out"]:
        out = odict["out"][0]
        if out[-1] != "/":
            out += "/"
        if out[0] == "/":
            out = out[1:]
        odict["out"] = out    
    # Returns the input/output arguments as a labeled array
    return odict


def get_elapsed_time(start=0):
    end = time.time() - start
    time_label = 'seconds'    
    if end > 3600:
        end = end/3600
        time_label = 'hours'
    elif end > 60:
        end = end/60
        time_label = 'minutes'
    return end, time_label


#v This function locates the relevant input data by input directory and file extension:
def find_input_data(indir, suffix, zeros_only=True, sort_time=True):
    if zeros_only:                                                  #< 'zeros_only' is currently hard-coded as we only want to work with the
        file_list = glob.glob(indir + "*0000" + suffix)             #< highest resolution filterbank files, which are labeled as '*.0000.fil'
    else:
        file_list = glob.glob(indir + "*" + suffix)                 #< This line does nothing at present
    assert len(file_list) >= 1                                      #< This is a check to make sure files were indeed found and listed
    if sort_time:                                                     #< Currently hard-coded so that our list of files is ordered by file_id
        file_list = sort_by_timestamp(file_list)
    return file_list                                                #< Returns the list of filterbank .fil files in the input directory to work on


def cull_dats(datfile,min_drift):                                                    
    file_array = open(datfile, 'r').read().splitlines()
    head_num = len(file_array) - 9            
    header = pd.read_csv(datfile,skipfooter=head_num, engine='python',index_col=False)
    headlist = header.values.tolist()
    df=pd.read_csv(datfile, sep=r'\t', header=[7],skiprows=[8], engine='python')
    result = df[abs(df['Drift_Rate ']) > min_drift].reset_index(drop=True)
    with open(datfile,'w') as ict:
        dummyvar = ict.write(str(headlist[1])[2:-2])
        dummyvar = ict.write('\n')
        for line in headlist:
            dummyvar = ict.write('\t'.join(line))
            dummyvar = ict.write('\n')
    with open(datfile, 'a') as f:
        for i in range(0,len(result)):
            for k in range(0,len(result.columns)):
                j = result.columns[k]
                if ((j == '# Top_Hit_# ') or (k == 'Coarse_Channel_Number ') or (k == 'Full_number_of_hits') or (k=='SEFD ')):
                    dummyvar = f.write(result[j][i].astype(str) + '\t')
                else:
                    dummyvar = f.write((f'{result[j][i]:.6f}') + '\t')
            dummyvar = f.write('\n')
    print(f'\nDrift rates below {min_drift:.4f} have been culled from ' + datfile)
    return None


#v This function converts the .fil files to .h5 files with blimpy:
def convert_to_h5(file, outdir="./", clobber=False):
    pre, ext = os.path.splitext(os.path.basename(file))                 #< Split the input filename by its file extension
    in_file = os.path.dirname(file) + "/" + pre + ".h5"
    out_file = outdir + pre + ".h5"                                     #< Create the output filename with the new extension .h5
    if os.path.isfile(out_file):                                        #< Check to see if the output filename already exists
        if not clobber:                                                 #< If the file exists and there is no clobber command then notify
            print("\n" + pre + ".h5" + " already exists. Moving on...")     
            return out_file                                             #< and return the .h5 filename without converting the .fil again
        else:                                                           #< Otherwise...
            os.remove(out_file)                                         #< If clobber is true, remove the existing output file to overwrite  
            print("%s has been removed successfully" %out_file)                           #< Notification of the completed process
    elif os.path.isfile(in_file):                                       # If it's not in the outdir, but it is in the indir...
        if not clobber:
            # move h5 file to output folder temporarily
            mv_h5_out = 'mv ' + os.path.basename(in_file) + " " + outdir
            call([mv_h5_out],shell=True)    
            return out_file
        else:
            os.remove(in_file)
            print("%s has been removed successfully" %out_file)                           #< Notification of the completed process
    print("\nConverting " + pre + " to HDF5...")                           #< Notification that the process will commence
    print("\nPlease ignore the following warning and INFO loggers from blimpy...")
    bp.fil2h5.make_h5_file(file, out_dir=outdir)                        #< Pass the .fil file through blimpy to make the .h5 file
    print("\nFile complete:\n" + out_file)       
    return out_file                                                     #< Return the .h5 output file


#v This function runs turboSETI, searching the .h5 input file and populating a .dat file 
#v with hits based on the specified max drift rate and signal to noise ratio:
def run_turbo_seti(file, max_drift=np.nan, min_drift=0, min_snr=10.0, outdir="./", clobber=False, gpu_backend=False):
    assert max_drift != np.nan                                              #< Make sure the drift rate input is a number
    pre, ext = os.path.splitext(os.path.basename(file))                     #< Split the input filename by its file extension
    out_file = outdir + pre + ".dat"                                        #< Create the output filename with the new extension .dat
    if os.path.isfile(out_file):                                            
        if not clobber:
            print("\n" + pre + ".dat" + " already exists. Moving on...")
            return out_file
        else:
            print("Removing old .dat to make way for the new one...")
            os.remove(out_file)
            out_log = outdir + pre + ".log"
            os.remove(out_log)
            print("%s has been removed successfully\n" %out_file) 
    # call FindDoppler
    log_level_int=logging.WARNING       # lessen the verbosity of the output
    fdop = FindDoppler(datafile=file, max_drift=max_drift, min_drift=min_drift, snr=min_snr, 
                       out_dir=outdir, gpu_backend=gpu_backend, log_level_int=log_level_int)
    # search for hits and report elapsed time.
    startt = time.time()
    fdop.search()
    delt, time_label = get_elapsed_time(startt)
    print(f"\nfind_doppler.py populated the corresponding .dat file with hits in %.2f {time_label}." %delt)
    print("\n")
    # return the .dat file name
    return out_file
    

def main():
    print("\nExecuting program...")
    start = time.time()

    # parse any command line arguments
    cmd_args = parse_args()
    indir = cmd_args["indir"]
    clobber = cmd_args["clobber"]
    max_drift_factor = cmd_args["MaxDF"]
    min_drift_factor = cmd_args["MinDF"]
    min_SNR = cmd_args["SNR1"]
    gpu_backend = cmd_args["GPU"]
    outdir = cmd_args["out"]

    # Make sure the data directory is set correctly or else exit program
    if not os.path.isdir(indir):
        print("\n Specified directory does not exist. Exiting... \n")
        sys.exit()

    # deal with GPU stuff and set output directory either way
    if gpu_backend:
        import cupy
        outdir = indir + outdir + "_gpu/"
    else:
        outdir = indir + outdir

    # make the output directory if it doesn't already exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # get appropriate list of .fil files in directory to work on
    fil_list = find_input_data(indir, '.fil')
    lst_fil = outdir + "fil_files.lst"
    # purge any existing fil_files.lst before creating the new one
    if os.path.isfile(lst_fil):
        os.remove(lst_fil)
        print("%s has been detected and purged to make way for the new list of '.fil' files." %lst_fil) 
    with open(lst_fil, 'w') as f:
        for item in fil_list:
            f.write("%s\n" % item)
    assert os.path.isfile(lst_fil)

    # loop over input files, reduce to .h5 and comb for hits into .dat file
    print("\nPlease wait ...")
    drifts=[]
    fil_array=[]
    start_dop = time.time()
    for fil in fil_list:
        # convert filterbank file to HDF5
        print(f"\nWorking on file number %s in the list..." %(fil_list.index(fil)+1))
        h5file = convert_to_h5(fil, outdir=outdir) 
        # Get max/min drift rates based on target barycentric drift and user input drift factors
        print("\nGetting max drift rate based on target coordinates...\n")
        drift_rate = get_drift_rate(h5file)
        drifts.append(drift_rate)
        fil_array.append(get_file_id(fil))
        print(f"\nBarycentric drift rate (pre factor) is {drift_rate:.2f}")
        max_drift_rate = drift_rate * max_drift_factor
        min_drift_rate = drift_rate * min_drift_factor                                              
        print(f"Using a drift rate range between {min_drift_rate:.2f} and {max_drift_rate:.2f}.")   
        # call FindDoppler
        datfile = run_turbo_seti(h5file,
                                 min_snr=min_SNR,
                                 outdir=outdir,
                                 clobber=clobber,
                                 max_drift=max_drift_rate,
                                 gpu_backend=gpu_backend)
        # move h5 file back to input folder
        mv_h5_in = 'mv ' + outdir + os.path.basename(h5file) + " " + indir
        call([mv_h5_in],shell=True)
        # apply minimum drift rate filter to dat files if indicated
        if min_drift_rate:
            cull_dats(datfile=datfile,min_drift=min_drift_rate)
            print('\nCulling complete.')
        else:
            min_drift_factor=0

    #This block prints the elapsed time for the entire find_doppler sequence.
    end_dop, time_label_dop = get_elapsed_time(start_dop)
    print(f"\nTotal time to populate dats with hits using find_doppler: %.2f {time_label_dop}.\n" %end_dop)
    
    # This block populates a 'dat_files.lst' file with a list of the .dat files in the output directory
    # First, remove any existing dat list    
    print("\nRemoving old dat lists...\n")
    old_files = glob.glob(outdir + "dat*.lst")
    if not old_files:
        print("There are no old .lst files to purge. Moving on...")
    else:
        for file in old_files:
            os.remove(file)
            print("%s has been removed successfully" %file) 
    # Next, make the new dat list
    dat_list = find_input_data(outdir, '.dat')
    lst_dat = outdir + "dat_files.lst"
    if os.path.isfile(lst_dat):
        os.remove(lst_dat)
    with open(lst_dat, 'w') as f:
        for item in dat_list:
            f.write("%s\n" % item)

    # This block calculates and prints the max/min drift rate range.
    max_max_drift_rate=max(drifts) * max_drift_factor
    min_min_drift_rate=min(drifts) * min_drift_factor
    print(f'\nAll hits within drift rate range of {min_min_drift_rate:.4f} and {max_max_drift_rate:.4f} Hz/s.')
    
    # This block prints an array of all the barycentric drift rates for each filterbank file.    
    formatted_drifts = [ '%.4f' % elem for elem in drifts]  
    print(f"\nBarycentric drift rates (pre factor) are {formatted_drifts}")
    print(f"\nThese correspond to the files with file ids: {fil_array}")

    # This block prints the elapsed time of the entire program.
    end, time_label = get_elapsed_time(start)
    print(f"\nTotal time to execute this program: %.2f {time_label}.\n" %end)
    return None
# run it!
if __name__ == "__main__":
    main()
