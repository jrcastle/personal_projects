import os


f = open("files3.txt")

fileList = f.readlines()
Nfiles = len(fileList)
count = 1

for line in fileList:
    line = line.replace("\n","")
    line = line.replace("/lio/lfs/cms","")
    fName = line.split("/")[9]
    if os.path.exists(fName):
        command= "rm " + str(fName)
        os.system(command)
    command = "xrdcp root://cmsxrootd.fnal.gov/" + str(line) + " " + str(fName)
    print "Transferring file " + str(fName) + " ...\t" + str(100.*float(count)/Nfiles) + "% Complete"
    os.system(command)
    count = count + 1
#END FOR

print "DONE!"
