##-----------------------------------------------------------------------------------------##
##-- A script that uses the DAS client to split a PD's list of files into separate lists --##
##-- based on pre-defined IOVs. das_client.py must be in the working directroy!          --##
##-----------------------------------------------------------------------------------------##
##-- Author: James Castle 

import os

##- Define PD, IOVs, dummy file, and output files.
PD     = "/HIHardProbes/HIRun2015-TkAlMinBiasHI-PromptReco-v1/ALCARECO"
IOV1   = 263203
fDummy = "RunNumbers.txt"
f1     = str(PD.split("/")[0]) + "_" + "Runs1-" + str(IOV1 - 1) + "_" + str(PD.split("/")[-1]) + ".py"
f2     = str(PD.split("/")[0]) + "_" + "Runs" + str(IOV1) + "-Inf_" + str(PD.split("/")[-1]) + ".py" 

##-- Check to see if output and dummy files exist in the working directory.  If so, delete them.
if os.path.exists(fDummy):
    command = "rm " + str(fDummy)
    print command
    os.system(command)

if os.path.exists(f1):
    command = "rm " + str(f1)
    print command
    os.system(command)

if os.path.exists(f2):
    command = "rm " + str(f2)
    print command
    os.system(command)

##
##-- IOV1
##

##-- Get all runs for this IOV from DAS and pipe into a dummy file
command="./das_client.py --query=\"run dataset=" + str(PD) + " | grep run.run_number | sort run.run_number\" --limit=0 | awk \'$1<"+str(IOV1)+"\'> RunNumbers.txt"
print command
os.system(command)

##-- Loop over all runs in the IOV and pipe their files into a list that will be converted to a python cff
f = open("RunNumbers.txt")
for line in f.readlines():
    line=line.replace("\n","")
    command="./das_client.py --query=\"file run=" + str(line) + " dataset=" + str(PD) + "\" --limit=0 >> " + str(f1)
    os.system(command)
##--End For

f.close()
##-- Remove the dummy file so that the next IOV is not contaminated
command = "rm RunNumbers.txt"
os.system(command)

##
##-- IOV2
##

##-- Get all runs for this IOV from DAS and pipe into a dummy file
command="./das_client.py --query=\"run dataset=" + str(PD) + " | grep run.run_number | sort run.run_number\" --limit=0 | awk \'$1>="+str(run)+"\'> RunNumbers.txt"
print command
os.system(command)

##-- Loop over all runs in the IOV and pipe their files into a list that will be converted to a python cff 
f = open("RunNumbers.txt")
for line in f.readlines():
    line=line.replace("\n","")
    command="./das_client.py --query=\"file run=" + str(line) + " dataset=" + str(PD) + "\" --limit=0 >> " + str(f2)
    os.system(command)
##--End For

f.close() 

print "DONE!"
