#!/usr/bin/env python
#################################################################################################
# Dependencies
#################################################################################################
import os
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d',    '--dir',       help='path to UFS Regression Test output',                           required=True)
parser.add_argument('-n',    '--case_name', help='name of case',                                                 required=True)
parser.add_argument('-lonl', '--lon_range', help='longitude range, separated by a space', nargs=2,   type=int, required=True)
parser.add_argument('-latl', '--lat_range', help='latitude range, separated by a space',  nargs=2,   type=int, required=True)
parser.add_argument('-dlat', '--dlat',      help='latitude spacing',                                 type=int, required=True)
parser.add_argument('-dlon', '--dlon',      help='longitude spacing',                                type=int, required=True)
parser.add_argument('-dt',   '--timestep',  help='SCM timestep, in seconds',                         type=int,   required=True)
parser.add_argument('-cres', '--C_RES',     help='UFS spatial resolution',                           type=int,   required=True)
parser.add_argument('-sdf',  '--suite',     help='CCPP suite definition file to use for ensemble',               required=True)
parser.add_argument('-fout', '--fileOUT',   help='Text file containing lon/lat points',                        required=False, default='lonlats.txt')

def main():

    # Get command line arguments
    args  = parser.parse_args()
    dir = args.dir+'/'+args.case_name

    com_pre = './UFS_forcing_ensemble_generator.py -d '+dir+' -sc --C_RES '+str(args.C_RES)+' -dt '+str(args.timestep)+' -n '+args.case_name+' -sdf '+args.suite

    lons = []
    lats = []
    count = 0
    for lon in range(args.lon_range[0],args.lon_range[1],args.dlon):
        for lat in range(args.lat_range[0],args.lat_range[1],args.dlat):
            count = count + 1
            lons.append(lon)
            lats.append(lat)
        # end for
    # end for

    comA =  '-lons '
    for ij in lons:
        comA = comA + str(ij)+' '
    # end for
    comB = '-lats '
    for ij in lats:
        comB = comB + str(ij)+' '
    # end for
    com = com_pre+comA+comB
    print(com)
    print("#######################################################################################")

    # For python list input file
    comA = ''
    countA = 0
    for ij in lons:
        countA = countA + 1
        if (countA < count):
            comA = comA + str(ij)+', '
        else:
            comA = comA + str(ij)
        #end if
    # end for

    comB = ''
    countB = 0 
    for ij in lats:
        countB = countB + 1
        if (countB < count):
            comB = comB + str(ij)+', '
        else:
            comB = comB + str(ij)
        # end if
    # end for

    fileID = open(args.fileOUT,'w')
    fileID.write(comA+"\n")
    fileID.write(comB+"\n")
    fileID.close()
    com = com_pre+' -fxy '+args.fileOUT
    print(com)

if __name__ == '__main__':
    main()
