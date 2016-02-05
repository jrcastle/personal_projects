##---------------------------------------------------------------------------------------------------##
##-- A simple script that reads a file that contains a list of all files from a PD taken from DAS  --##
##-- and randomly splits these files into two separate lists.  The format in which these are split --##
##-- is for Millepede jobs in CMS, but can easily be adapted to fit the user's needs               --##
##---------------------------------------------------------------------------------------------------##
##-- Author: James Castle

import os
import random

random.seed()

##-- Output files
file1 = "HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"
file2 = "HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"

##-- Check if these files are within the working directory currently.  If so, remove them
if os.path.exists(file1):
    command = "rm " + str(file1)
    print command
    os.system(command)

if os.path.exists(file2):
    command = "rm " + str(file2)
    print command
    os.system(command)

##-- Open the master list
file = open("temp.txt")

#Add necessary CastorPool for Millepede data lists
command = "echo \"CastorPool=cmscaf\" >> " + str(file1)
print command
os.system(command)
command = "echo \"CastorPool=cmscaf\" >> " + str(file2)
print command
os.system(command)

##-- Loop over the files in the master list and throw a random number for each file.
##-- Use the RNG to distribute files as you please
for line in file.readlines():
    line=line.replace("\n","")

    i = random.randrange(1,100,1)

    if i <= 7:
        command = "echo \"" + str(line) + "\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_ALIGN.txt"
    else:
        command = "echo \"" + str(line) + "\" >> HIHardProbes_HIRun2015-TkAlMinBiasHI-PromptReco-v1_ALCARECO_VALIDATE.txt"

    os.system(command)
##-- END FOR

print "DONE!"
