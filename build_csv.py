# This program combs through the .dat files and filters the hits based on the SNR and filter_threshold inputs 
# and outputs the resulting hits in .csv file.

import numpy as np
import os
import argparse
import glob
import time
import turbo_seti
import blimpy as bp
from parse_filname import sort_by_timestamp
from subprocess import call
from turbo_seti.find_event.find_event_pipeline import find_event_pipeline

def parse_args():
    parser = argparse.ArgumentParser(description='Process BL filterbank data.')
    parser.add_argument('fildir', metavar='/fil_file_directory/', type=str, nargs=1,
                        help='directory containing the .fil and .h5 files')
    parser.add_argument('datf', metavar='/dat_file_subfolder/', type=str, nargs=1,
                        help='subfolder within fildir containing the .dat files')
    parser.add_argument('-GPU', action='store_true',
                        help='use GPU acceleration if possible')
    parser.add_argument('-s2', dest='SNR2', metavar='#Number', type=float, default=10,
                        help='Set the signal to noise ratio as a number to search for hits')
    parser.add_argument('-f', dest='filter_threshold', metavar='1, 2, or 3', type=int, default=1, choices=range(1, 4),
                        help='Set the filter threshold to 1, 2, or 3')
    parser.add_argument('-n', dest='number_in_cadence', metavar='#total', type=int, default=6,
                        help='The number of files in each cadence. Default set to max of 6.')
    parser.add_argument('-t', metavar='ON source name', type=str, nargs=1,
                        help='ON source string in filename')
    args = parser.parse_args()
    # Check for trailing slash in the directory path and add it if absent
    odict = vars(args)
    if odict["fildir"]:
        indir = odict["fildir"][0]
        if indir[-1] != "/":
            indir += "/"
        odict["fildir"] = indir  
    if odict["datf"]:
        indir = odict["datf"][0]
        if indir[-1] != "/":
            indir += "/"
        if indir[0] == "/":
            indir = indir[1:]
        odict["datf"] = indir  
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


def run_find_event_pipeline(fildir, datdir, datlist, SNR, filter_threshold=np.nan, number_in_cadence=np.nan,
                                       on_source_complex_cadence=False):
    assert filter_threshold != np.nan
    assert number_in_cadence != np.nan
    # write to .lst, as required by find_event_pipeline
    lst_dat = datdir + "dat_files.lst"
    if not os.path.isfile(lst_dat):
        with open(lst_dat, 'w') as f:
            for item in datlist:
                f.write("%s\n" % item)
    # move the h5 files to the dat directory
    for file in datlist:
        h5_file = fildir + os.path.basename(file).replace('.dat','.h5')
        if not os.path.isfile(h5_file):
            fil = os.path.basename(h5_file).replace('.h5','.fil')
            bp.fil2h5.make_h5_file(fil, out_dir=datdir)
        else:
            mvh5_in = 'mv ' + h5_file + " " + datdir
            call([mvh5_in],shell=True)
    # run find_event_pipeline
    print("\nRunning find event pipeline...")
    find_event_pipeline(lst_dat,
                        SNR_cut=SNR, 
                        check_zero_drift=False, 
                        filter_threshold=filter_threshold,
                        number_in_cadence=number_in_cadence,
                        on_source_complex_cadence=on_source_complex_cadence)
    # move the h5 files back to the fil directory to save space
    for file in datlist:
        mvh5_out = 'mv ' + datdir + os.path.basename(file).replace('.dat','.h5') + " " + fildir
        call([mvh5_out],shell=True)


def main():
    print("\nExecuting program...")
    start = time.time()

    # parse any command line arguments
    cmd_args = parse_args()
    fildir = cmd_args["fildir"]
    datf = cmd_args["datf"]
    SNR = cmd_args["SNR2"]
    gpu_backend = cmd_args["GPU"]
    filter_threshold = cmd_args["filter_threshold"]
    number_in_cadence = cmd_args["number_in_cadence"]
    target_name = str(cmd_args["t"])[2:-2]
    datdir=fildir+datf
    dat_list = find_input_data(datdir, '.dat')
    if target_name:
        number_in_cadence = len(dat_list)
    else:
        target_name=bool(0)
	# create the csv file
    run_find_event_pipeline(fildir=fildir,
                            datdir=datdir,
                            datlist=dat_list,
                            SNR=SNR,
                            filter_threshold=filter_threshold,
                            number_in_cadence=number_in_cadence,
                            on_source_complex_cadence=target_name)
    # print elapsed time
    end, time_label = get_elapsed_time(start)
    print(f"\nTotal time to execute this program: %.2f {time_label}.\n" %end)
    return None
# run it!
if __name__ == "__main__":
    main()