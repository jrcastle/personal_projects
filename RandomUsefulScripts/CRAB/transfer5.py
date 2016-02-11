##------------------------------------------------------------------------------------------##
##-- Follow-up script to splitPD_Vandy.py.  This script reads in a single split file,     --##
##-- generates an xrdcp command specific to an element in the list, and then executes the --##
##-- command. The user needs to have a cmsenv already set and a proxy in place.           --##
##------------------------------------------------------------------------------------------##
##-- Author: James Castle
import os

#-- Open one of the split files
f = open("files5.txt")

##-- Count the number of files to be transferred to track progress
fileList = f.readlines()
Nfiles = len(fileList)
count = 1

##-- Loop over the list of files, generate the xrdcp command, and execute it
for line in fileList:
    line = line.replace("\n","")
    line = line.replace("/lio/lfs/cms","")
    fName = line.split("/")[9]
    ##-- Check to see if the file is already in the working directory.  If so, delete it.  Else xrdcp will complain
    if os.path.exists(fName):
        command = "rm " + str(fName)
        os.system(command)
    command = "xrdcp root://cmsxrootd.fnal.gov/" + str(line) + " " + str(fName)
    print "Transferring file " + str(fName) + " ...\t" + str(100.*float(count)/Nfiles) + "% Complete"
    os.system(command)
    count = count + 1
#END FOR

print "DONE!"
