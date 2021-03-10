# This code will take a typical set of filterbank data observed in an on/off cadence
#       and run it through turboSETI and blimpy to produce waterfall plots.
#
# This requires that:
#       the files be labeled with the standard BL at GBT filename convention.
#       the observations are ordered in their on/off cadence.
#       the filterbank files are spliced together.
#       the folder containing the filterbank files contain only the .fil files
#               pertaining to the on/off observation and no calibration files.
#               All other .fil files should be moved to a different folder.
#
# This program runs through the following steps:
#       1. Obtains coordinates from the headers of the filterbank files using blimpy.
#       2. Calculates the maximum drift rate to search for based on barycentric motion
#               multiplied by an optional input factor.
#       3. Uses the blimpy package to convert the .fil files to .h5 files.
#       4. Runs find_doppler in turbo_seti to comb through the h5 files and populate a
#               list of hits in a .dat file based on the max drift rate, and a
#               signal to noise ratio as an optional input, which defaults to 10.
#       5. Runs find_event_pipeline in turbo_seti to further filter the hits listed in
#               the .dat files into a .csv file based on another optional signal to
#               noise input value that defaults to 10, and a filter_threshold value
#               that is explained below.
#       6. Runs plot_event_pipeline to cross referene the .csv file and the original
#               filterbank data files and produce waterfall plots of any signals.
#
# The filter_threshold is meant to serve as automatic RFI rejection and can have one
#       of three input values:
#       1. Hits are filtered based solely on SNR.
#       2. Hits are filtered on SNR and must occur in at least one 'ON' source and no
#               'OFF' sources, as dictated the observing cadence.
#       3. Hits are filtered on SNR and must occur in every 'ON' source as well as
#               no 'OFF' sources, as dictated by the observing cadence.
#
# Note that level 1 provides little to no RFI rejection and level 3 is the most strict,
#       requiring a positive detection to appear multiple times. Level 2 is recommended
#       for new users.
#
# This code will output .h5, .dat, and .log files into an output '/processed' folder,
#       as well as a list of .fil files used in a fil_files.lst file. It will also
#       output a list of the .dat files in a dat_files.lst, a .csv file with the refined
#       table of hits, and any waterfall plots it makes into another subfolder.
#
# The input arguments for the code are currently:
#       '/input_directory/' where all the filterbank .fil files are located, with or
#               without the trailing slash. Note: this is the only mandatory argument.
#       '-clobber' if present this will purge the output directory of all files in
#               order to overwrite them, or skip past any that exist if it is left off.
#       '-GPU' utilization is technically possible with turboSETI, though currently
#               untested by myself. I intend to test and update this later. Default off.
#       '-d' this sets the factor by which you want to multiply the max drift rate.
#               The default value for this is 1.
#       '-s1' this sets the first minimum signal to noise level used by find_doppler.py
#               when scanning the data for hits. Default value is set to 10.
#       '-s2' this sets the second signal to noise level used by find_event_pipeline to
#               further filter the list of hits. Default value is set to 10.
#       'n' this sets the total number of files in the observation cadence, including
#               both on and off files. The maximum value for this is the default at 6.
#       'f' this sets the filter_threshold as described above. Default is 1.
#
# As an example, this code can be run on a pre-curated folder of .fil files from the
#       the terminal window like so:
# > python nodding_cadence_pipeline.py /input/path/ -clobber -m 10 -s2 20 -n 6 -f 2
#

#v PACKAGES:
import numpy as np
import os
import sys
import pdb
import glob
import time
import argparse
import csv
import logging
import pandas as pd
import turbo_seti
import blimpy as bp
from subprocess import call
from seti_lens import get_drift_rate
# from input_Tabby import indir, clobber, MDR, SNR
from turbo_seti.find_doppler.find_doppler import FindDoppler
from turbo_seti.find_event.find_event_pipeline import find_event_pipeline
from turbo_seti.find_event.plot_event_pipeline import plot_event_pipeline

#v The pickle package and function are for code testing purposes only:
import pickle
def my_little_pickle(thingy, name):
    print(f"\nMy little pickle: {name}\n{pickle.loads(pickle.dumps(thingy))}\n")

#v FUNCTIONS:
#v This function parses the input arguments:
def parse_args():
    parser = argparse.ArgumentParser(description='Process BL filterbank data.')
    parser.add_argument('-p', dest='par', metavar='input_parameters_filename', type=str, nargs=1,
                        help='filename containing the input variables without extension')
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
    parser.add_argument('-s2', dest='SNR2', metavar='#Number', type=float, default=10,
                        help='Set the signal to noise ratio as a number to search for hits')
    parser.add_argument('-f', dest='filter_threshold', metavar='1, 2, or 3', type=int, default=1, choices=range(1, 4),
                        help='Set the filter threshold to 1, 2, or 3')
    parser.add_argument('-n', dest='number_in_cadence', metavar='#total', type=int, default=6,
                        help='The number of files in each cadence. Default set to max of 6.')
    parser.add_argument('-out', metavar='/outut_directory/', type=str, default='processed/', nargs=1,
                        help='directory to put output files into')
    parser.add_argument('-t', metavar='ON source name', type=str, nargs=1,
                        help='ON source string in filename')
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
def convert_to_h5(file, indir="./", clobber=False):
    pre, ext = os.path.splitext(os.path.basename(file))                 #< Split the input filename by its file extension
    out_file = indir + pre + ".h5"                                     #< Create the output filename with the new extension .h5
    if os.path.isfile(out_file):                                        #< Check to see if the output filename already exists
        if not clobber:                                                 #< If the file exists and there is no clobber command then notify
            print("\n" + pre + ".h5" + " already exists. Moving on...")
            return out_file                                             #< and return the .h5 filename without converting the .fil again
        else:                                                           #< Otherwise...
            os.remove(out_file)                                         #< If clobber is true, remove the existing output file to overwrite
            print("%s has been removed successfully" %out_file)
    print("\nConverting " + pre + " to HDF5...")                           #< Notification that the process will commence
    print("\nPlease ignore the following warning and INFO loggers from blimpy...")
    bp.fil2h5.make_h5_file(file, out_dir=indir)                        #< Pass the .fil file through blimpy to make the .h5 file
    print("\nFile complete:\n" + out_file)                              #< Notification of the completed process
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
    log_level_int=logging.WARNING
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


def run_find_event_pipeline(datdir, SNR, filter_threshold=np.nan, number_in_cadence=np.nan,
                            on_source_complex_cadence=False):
    assert filter_threshold != np.nan
    assert number_in_cadence != np.nan
    #get list of files and sort by file ID
    dat_list = glob.glob(datdir + "*.dat")
    dat_list = sort_by_timestamp(dat_list)
    # write to .lst, as required by find_event_pipeline
    lst_dat = datdir + "dat_files.lst"
    if os.path.isfile(lst_dat):
        os.remove(lst_dat)
    with open(lst_dat, 'w') as f:
        for item in dat_list:
            f.write("%s\n" % item)
    assert os.path.isfile(lst_dat)

    # construct csv file name
    csv_name = datdir + "hits_SNR_" + SNR + "_f_" + str(filter_threshold) + ".csv"
    if os.path.isfile(csv_name):
        os.remove(csv_name)

    # run find_event_pipeline
    print("\nRunning find event pipeline...")
    find_event_pipeline(lst_dat,
                        SNR_cut=SNR,
                        csv_name=csv_name,
                        check_zero_drift=False,
                        filter_threshold=filter_threshold,
                        number_in_cadence=number_in_cadence,
                        on_source_complex_cadence=on_source_complex_cadence)
    # try:
    #     csv_name = (glob.glob('*f' + str(filter_threshold) + '_snr' + str(SNR) + '*.csv'))[0]
    # except:
    #     csv_name = []
    # print("\n")
    return csv_name


def main():
    print("\nExecuting program...")
    start = time.time()
    # parse any command line arguments
    cmd_args = parse_args()
    if not cmd_args["par"]:
        # get the input variables from the command arguments if no input variables file listed
        indir = cmd_args["indir"]
        clobber = cmd_args["clobber"]
        max_drift_factor = cmd_args["MaxDF"]
        min_drift_factor = cmd_args["MinDF"]
        min_SNR = cmd_args["SNR1"]
        gpu_backend = cmd_args["GPU"]
        filter_threshold = cmd_args["filter_threshold"]
        number_in_cadence = cmd_args["number_in_cadence"]
        SNR = cmd_args["SNR2"]
        outdir = cmd_args["out"]
        target_name = str(cmd_args["t"])[2:-2]
    else:                               # still trying to get the input file working
        parameters = (cmd_args["par"])
        my_little_pickle(parameters)
        v = globals().update(importlib.import_module(parameters).__dict__)
        my_little_pickle(v)
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
    # make the "processed" directory if it doesn't already exist
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # get appropriate list of .fil files in directory to work on
    fil_list = find_input_data(indir, '.fil')
    # write .fil files to .lst file for plot_event_pipeline later
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
        h5_file = convert_to_h5(fil, indir=indir)

        # Get max drift rate
        print("\nGetting max drift rate based on target coordinates...\n")
        drift_rate = get_drift_rate(h5_file)
        drifts.append(drift_rate)
        fil_array.append(get_file_id(fil))
        print(f"\nBarycentric drift rate (pre factor) is {drift_rate:.2f}")
        max_drift_rate = drift_rate * max_drift_factor
        min_drift_rate = drift_rate * min_drift_factor
        print(f"Using a drift rate range between {min_drift_rate:.2f} and {max_drift_rate:.2f}.")

        # run turbo_seti
        datfile = run_turbo_seti(h5_file,
                                 min_snr=min_SNR,
                                 outdir=outdir,
                                 clobber=clobber,
                                 max_drift=max_drift_rate,
                                 gpu_backend=gpu_backend,
                                 clobber=clobber)

        # workaround for broken min drift rate
        if min_drift_rate:
            cull_dats(datfile=datfile,min_drift=min_drift_rate)
            print('\nCulling complete.')
        else:
            min_drift_factor=0

        pdb.set_trace()

    # multiply drift rates by factor and print
    max_max_drift_rate=max(drifts) * max_drift_factor
    min_min_drift_rate=min(drifts) * min_drift_factor
    print(f'\nAll hits within drift rate range of {min_min_drift_rate:.4f} and {max_max_drift_rate:.4f} Hz/s.')

    formatted_drifts = [ '%.4f' % elem for elem in drifts]
    print(f"\nBarycentric drift rates (pre factor) are {formatted_drifts}")
    print(f"\nThese correspond to the files with file ids: {fil_array}")

    end_dop, time_label_dop = get_elapsed_time(start_dop)
    print(f"\nTotal time to populate dats with hits using find_doppler: %.2f {time_label_dop}.\n" %end_dop)

    # remove old dat lists even if clobber off
    print("\nRemoving old dat lists...\n")
    old_files = glob.glob(outdir + "dat*.lst")
    if not old_files:
        print("There are no old .lst files to purge. Moving on...")
    else:
        for file in old_files:
            os.remove(file)
            print("%s has been removed successfully" %file)
    # get appropriate list of .dat files in directory
    dat_list = find_input_data(outdir, '.dat')

    # call run_find_event_pipeline and make csv file to be placed in the processed folder
    # os.chdir(outdir)               # change into output dir

    if target_name:
        number_in_cadence = len(dat_list)
    else:
        target_name=bool(0)

    csv_file = run_find_event_pipeline(datdir=outdir, SNR=SNR,
                                       filter_threshold=filter_threshold,
                                       number_in_cadence=number_in_cadence,
                                       on_source_complex_cadence=target_name)

    # now do the plotting
    # cross-reference the csv file above with the fil_files.lst
    if csv_file:
        # get the directory to save the plot to
        plot_dir = outdir + "f" + str(filter_threshold) + '_s' + str(int(SNR)) + '_plots/'
        if not os.path.isdir(plot_dir):
            os.mkdir(plot_dir)

        print("\nRunning plot_event_pipeline...")
        plot_event_pipeline(csv_file, lst_fil, plot_dir=plot_dir)

        # command = 'mv ' + indir + 'f' + str(filter_threshold) + '*.p?? ' + plot_dir
        # call([command],shell=True)

    print(f"\nBarycentric drift rates (pre factor) are {formatted_drifts}")
    print(f"\nThese correspond to the files with file ids: {fil_array}")

    end, time_label = get_elapsed_time(start)
    print(f"\nTotal time to execute this program: %.2f {time_label}.\n" %end)
    return None
# run it!
if __name__ == "__main__":
    main()
