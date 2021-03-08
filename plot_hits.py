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
from parse_filname import sort_by_timestamp
from turbo_seti.find_event.plot_event_pipeline import plot_event_pipeline

#v The pickle package and function are for code testing purposes only:
import pickle
def my_little_pickle(thingy, name):
    print(f"\nMy little pickle: {name}\n{pickle.loads(pickle.dumps(thingy))}\n") 

#v FUNCTIONS:
#v This function parses the input arguments:
def parse_args():
    parser = argparse.ArgumentParser(description='Process BL filterbank data.')
    parser.add_argument('fildir', metavar='/fil_file_directory/', type=str, nargs=1,
                        help='directory containing the .fil files')
    parser.add_argument('csvf', metavar='/csv_subfolder/', type=str, nargs=1,
                        help='folder within fildir that contains csv files to be plotted')
    args = parser.parse_args()
    # Check for trailing slash in the directory path and add it if absent
    odict = vars(args)
    if odict["fildir"]:
        indir = odict["fildir"][0]
        if indir[-1] != "/":
            indir += "/"
        odict["fildir"] = indir    
    if odict["csvf"]:
        indir = odict["csvf"][0]
        if indir[-1] != "/":
            indir += "/"
        if indir[0] == "/":
            indir = indir[1:]
        odict["csvf"] = indir 
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
    if sort_time:                                                     #< Currently hard-coded so that our list of files is ordered by file_id
        file_list = sort_by_timestamp(file_list)
    return file_list                                                #< Returns the list of filterbank .fil files in the input directory to work on


def main():
    print("\nExecuting program...")
    start = time.time()
    # parse any command line arguments
    cmd_args = parse_args()
    fildir = cmd_args["fildir"]
    csvf = cmd_args["csvf"]
    csvpath = fildir+csvf
    os.chdir(csvpath)
    lst_fil = csvpath + "fil_files.lst"
    csv_list = sorted(glob.glob('*.csv')) 
    # move the h5 files to the dat directory
    h5_indir = find_input_data(fildir,'.h5')
    h5_datdir = find_input_data(csvpath,'.h5')
    if not h5_datdir:
        if h5_indir:
            for file in h5_indir:
                h5_file_path = fildir + os.path.basename(file)
                mvh5_in = 'mv ' + h5_file_path + " " + csvpath
                call([mvh5_in],shell=True)
        else:
            fil_list = find_input_data(fildir,'.fil')
            for fil in fil_list:
                bp.fil2h5.make_h5_file(fil, out_dir=csvpath)

    if csv_list:
        for item in csv_list:
            csv_file = item         # csv string in plot_event_pipeline has to be filename only, not entire filepath.
            print("\nRunning plot_event_pipeline...")
            plot_event_pipeline(csv_file, 
                                lst_fil)
            plot_dir = csvpath + os.path.splitext(item)[0].split('.')[0] + '_plots/'
            if not os.path.isdir(plot_dir):
                os.mkdir(plot_dir)
            mv_plots = 'mv ' + fildir + os.path.splitext(item)[0].split('.')[0].split('_')[2] + '*.p?? ' + plot_dir  
            call([mv_plots],shell=True)
    else:
        print(f'\nThere is no csv file file in {csvpath}')

    # move the h5 files back to the fil directory to save space
    h5_done = find_input_data(csvpath,'.h5')
    for h5 in h5_done:
        mvh5_out = 'mv ' + csvpath + os.path.basename(h5) + " " + fildir
        call([mvh5_out],shell=True)

    end, time_label = get_elapsed_time(start)
    print(f"\nTotal time to execute this program: %.2f {time_label}.\n" %end)
    return None
# run it!
if __name__ == "__main__":
    main()