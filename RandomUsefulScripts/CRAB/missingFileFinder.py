import os

fOut = "MISSING.txt"
if os.path.exists(fOut):
    command = "rm " + fOut
    print command
    os.system(command)


fMaster = open("MASTER.txt")
mL = fMaster.readlines()
fCurrent = open("CURRENT.txt")
cL = fCurrent.readlines()

p = str(mL[0])
p = p.replace("\n", "")
sp = p.split("/")
path = ""
for i in range(0, len(sp)-1):
    path = path + str(sp[i]) + "/"

masterList = []
for line in mL:
    line = line.replace("\n", "")
    fName = line.split("/")[12]
    masterList.append(fName)

currentList = []
for line in cL:
    line = line.replace("\n", "")
    currentList.append(line)

for f in masterList:
    if f not in currentList:
        command = "echo \"" + path + f + "\" >> " + fOut
        print command
        os.system(command)
