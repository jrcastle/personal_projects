import os

if os.path.exists("RunNumbers.txt"):
    command = "rm RunNumbers.txt"
    print command
    os.system(command)

if os.path.exists("HIHardProbesALCARECO_Runs1-263202_VALIDATE_cff.py"):
    command = "rm HIHardProbesALCARECO_Runs1-263202_VALIDATE_cff.py"
    print command
    os.system(command)

if os.path.exists("HIHardProbesALCARECO_Runs263203-Inf_VALIDATE_cff.py"):
    command = "rm HIHardProbesALCARECO_Runs263203-Inf_VALIDATE_cff.py"
    print command
    os.system(command)

##-- IOV 1 - 263202
run = 263203
command="./das_client.py --query=\"run dataset=/HIHardProbes/HIRun2015-TkAlMinBiasHI-PromptReco-v1/ALCARECO | grep run.run_number | sort run.run_number\" --limit=0 | awk \'$1<"+str(run)+"\'> RunNumbers.txt"
print command
os.system(command)

f = open("RunNumbers.txt")
for line in f.readlines():
    line=line.replace("\n","")
    command="./das_client.py --query=\"file run=" + str(line) + " dataset=/HIHardProbes/HIRun2015-TkAlMinBiasHI-PromptReco-v1/ALCARECO\" --limit=0 >> HIHardProbesALCARECO_Runs1-263202_VALIDATE_cff.py"
    os.system(command)
##--End For

f.close()
command = "rm RunNumbers.txt"
os.system(command)

##-- IOV 263203 -Inf
command="./das_client.py --query=\"run dataset=/HIHardProbes/HIRun2015-TkAlMinBiasHI-PromptReco-v1/ALCARECO | grep run.run_number | sort run.run_number\" --limit=0 | awk \'$1>="+str(run)+"\'> RunNumbers.txt"
print command
os.system(command)

f = open("RunNumbers.txt")
for line in f.readlines():
    line=line.replace("\n","")
    command="./das_client.py --query=\"file run=" + str(line) + " dataset=/HIHardProbes/HIRun2015-TkAlMinBiasHI-PromptReco-v1/ALCARECO\" --limit=0 >> HIHardProbesALCARECO_Runs263203-Inf_VALIDATE_cff.py"
    os.system(command)
##--End For

f.close() 
