##----------------------------------------------------------------------------------------------------##
##-- A very simple python script that reads in a master list of files stored at Vanderbilt and      --##
##-- uses the lcg-del command to remove them.  This script is in it's first revision and it is      --## 
##-- planned to expand this to other storage elements.  This ONLY works for files, not directories! --##
##----------------------------------------------------------------------------------------------------##
##-- Author: James Castle
##-- Example to generate MASTER.txt:
##-- lcg-ls -b -D srmv2 'srm://se1.accre.vanderbilt.edu:6288/srm/v2/server?SFN=/lio/lfs/cms/store/user/<path-to-CRAB-output>' >> MASTER.txt

##-- IMPORTANT!!!!!!! DO NOT USE lcg-del ON FILES STORED AT MIT.  THEY WILL GET VERY MAD AT YOU!!!!!!
##-- THE PROPER PROTOCOL IS TO E-MAIL THE MIT CONTACT AND REQUEST FILES TO BE DELETED 

import os

##-- Open the master list
f = open("MASTER.txt")
files = f.readlines()


##-- Loop over the lines in the master file, generate delete command, and execute it
for line in files:
    line = line.replace("\n","")
    command = "lcg-del -D srmv2 -b -l srm://se1.accre.vanderbilt.edu:6288/srm/v2/server?SFN=" + str(line)
    print command
    os.system(command)
##-- END FOR

print "DONE!"







