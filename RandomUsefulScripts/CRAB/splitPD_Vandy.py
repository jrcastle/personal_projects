##-------------------------------------------------------------------------------------------------##
##-- A simple python script that reads in a master list of the output of an lcg-ls command from  --##
##-- Vanderbilt's T2.  The output contains all files from a CRAB job that have been staged out   --##
##-- and need to be retreived locally.  This script will then split the master into 5 lists that --##
##-- will be read by the transfer*.py scripts contained in this repository.  The end game is to  --##
##-- "parallelize" the transfer of files to the local machine and save time.  Once this script   --##
##-- is run, please proceed to run the transfer*.py scripts.                                     --##
##-------------------------------------------------------------------------------------------------##
##-- Author: James Castle
##--
##-- Sample lcg-ls command for Vanderbilt:
##-- lcg-ls -b -D srmv2 'srm://se1.accre.vanderbilt.edu:6288/srm/v2/server?SFN=/lio/lfs/cms/store/user/'

import os

##-- Output files
f1 = "files1.txt"
f2 = "files2.txt"
f3 = "files3.txt"
f4 = "files4.txt"
f5 = "files5.txt"

##-- Check if the output files exist in the working directory.  If so, remove them
if os.path.exists(f1):
    command = "rm " + str(f1)
    print command
    os.system(command)

if os.path.exists(f2):
    command = "rm " + str(f2)
    print command
    os.system(command)

if os.path.exists(f3):
    command = "rm " + str(f3)
    print command
    os.system(command)

if os.path.exists(f4):
    command = "rm " + str(f4)
    print command
    os.system(command)

if os.path.exists(f5):
    command = "rm " + str(f5)
    print command
    os.system(command)

##-- Open the master list
f = open("MASTER.txt")
fList = f.readlines()
NFiles = len(fList)
fCount = 1

##-- Loop through the master list and split the master list into 5 separate files
for line in fList:
    line=line.replace("\n","")

    fPercent = float( float( fCount ) / float( NFiles ) )

    if fPercent <= 0.2:
        command = "echo \"" + str(line) + "\" >> " + str(f1)
    if fPercent > 0.2 and fPercent <= 0.4:
        command = "echo \"" + str(line) + "\" >> " + str(f2)
    if fPercent > 0.4 and fPercent <= 0.6:
        command = "echo \"" + str(line) + "\" >> " + str(f3)
    if fPercent > 0.6 and fPercent <= 0.8:
        command = "echo \"" + str(line) + "\" >> " + str(f4)
    if fPercent > 0.8 and fPercent <= 1.0:
        command = "echo \"" + str(line) + "\" >> " + str(f5)

    os.system(command)
    fCount = fCount + 1
##-- ENDFOR

print "DONE!"
