#!/usr/bin/env python

#This is the main code that is run which queries DAS to get the latest run, finds the files in that run, starts the mille, pede jobs

import sys
import subprocess
import os
from ROOT import TH1F, TFile, TH1
from subprocess import *
import time
import shlex
import stat

#Tell the script whether it should be uploading results to the dropbox
uploadToDropbox = False #uploads them to 

skipExpressStreamFinishedCheck = True  #checks for run finished from tier0

sendEmailNotification = True

start = time.time()


print 'Start New Job ' + str(time.asctime()) 

#Get last run from file "ProcessedRuns.txt"
if os.path.exists("ProcessedRuns.txt"):
    ProcessedRunFile=open('ProcessedRuns.txt','r')
    ProcessedRuns = ProcessedRunFile.readlines()
    ProcessedRunFile.close()

    if len(ProcessedRuns) > 0:
        lastRun = ProcessedRuns[-1]
        lastRun = lastRun[:-1]
        print lastRun
    else:
        lastRun = "254986"
else:
    ProcessedRunFile=open('ProcessedRuns.txt','w')
    lastRun = "254986"
    
foundRun = False

if os.path.exists("RunNumbers.txt"): os.system("rm RunNumbers.txt")
    
command="./das_client.py --query=\"run dataset=/StreamHIExpress/HIRun2015-TkAlMinBiasHI-Express-v1/ALCARECO | grep run.run_number | sort run.run_number\" --limit=0 | awk \'$1>"+lastRun+"\'> RunNumbers.txt"
print command
os.system(command)
    
RunFile=open('RunNumbers.txt','r')
RunLines = RunFile.readlines()
RunFile.close()


if len(RunLines) > 0:
    for run_num in RunLines:
        command = "curl -k https://cmsweb.cern.ch/tier0/isexpresscomplete?run=" + run_num[:-1] + " >& curlOutput.txt"
        print command
        os.system(command)
        curlOut=open('curlOutput.txt','r')
        curlLine = curlOut.readlines()        
        isFinished = curlLine[-1]
        os.system("rm curlOutput.txt")


        
        runNumber = run_num[:-1]

        if isFinished == "True" or skipExpressStreamFinishedCheck:
            if os.path.exists("tempFiles.txt"):
                os.system("rm tempFiles.txt")
            callCommand = "./das_client.py --query=\"file run="+run_num[:-1]+" dataset=/StreamHIExpress/HIRun2015-TkAlMinBiasHI-Express-v1/ALCARECO | grep file.name\" --limit=0 > tempFiles.txt"
            print callCommand
            os.system(callCommand)
            
            # Text File with output of DAS query listing the file names
            FileNames=open('tempFiles.txt','r')
            ListFiles = FileNames.readlines()
            FileNames.close()
            #os.system("rm tempFiles.txt")
            
        
            if len(ListFiles) > 0:
                command = "./das_client.py --query=\"summary run="+run_num[:-1]+" dataset=/StreamHIExpress/HIRun2015-TkAlMinBiasHI-Express-v1/ALCARECO | grep summary.nevents\" --limit=0"
                p = subprocess.Popen(shlex.split(command),stdout = PIPE)
                output = p.stdout.read()
                Nevents = int(output.split('=')[-1][:-1])
                if Nevents < 10000: 
                    print "Fewer than 10,000 events in run, alignment not run for run", runNumber
                    print "Number of events: ", Nevents
		    continue
                command = "echo " + runNumber + " >> ProcessedRuns.txt"
                os.system(command)
                foundRun=True
                break
            else:
                print "No files for run number "+runNumber
        else:            
            print "Run Number "+runNumber+" is not finished with express stream"
            break
            

if not foundRun:
    sys.exit(0)

print runNumber

print "Number of files", len(ListFiles)
print "Number of events", Nevents

command = "./das_client.py --query=\"run="+runNumber+" | grep run.bfield\" --limit=0"
print command
p = subprocess.Popen(shlex.split(command),stdout = subprocess.PIPE)
bfield = p.stdout.read()
bfield = bfield.replace('\n','')

if float(bfield) < 3.7:
    print "WARNING 3.8T alignment used for Run " + runNumber + ", where bfield = "+str(bfield)
    command = "echo \"WARNING 3.8T alignment used for Run " + runNumber + ", where bfield = "+str(bfield)+"\" | mail -s \"WARNING 3.8T Alignment Used By Mistake\" jcastle@cern.ch"
    os.system(command)    

#Make directory for specific run and move files into there
directory = "Results" + runNumber

if os.path.exists(directory):
    command = "rm -rf "+directory
    os.system(command)

command = "mkdir " + directory
os.system(command)
command = "cp tempFiles.txt " + directory + "/."
os.system(command)
os.getcwd()
command = "cp *.sh " + directory + "/."
os.system(command)
os.getcwd()
command = "cp *.py " + directory + "/."
os.system(command)
os.getcwd()
os.chdir(directory)
os.getcwd()

#Run alignment
os.system('bash splitInput.sh')
#command_wait = "bsub -q 8nh -o output_pede.txt -e error_pede.txt -J Pede_2015 -w \""
command_wait = "bsub -q cmsexpress -o output_pede.txt -e error_pede.txt -J Pede_2015 -w \""
tempNum = min(80,len(ListFiles))
for file in range(0,tempNum):
	#print file
	if file != 0:
		command_wait = command_wait + " &&"
	#command = "bsub -q 8nh -o output_{1}.txt -e error_{1}.txt -J MinBias_2015_{1} automationMinBias.sh {0} {1}".format(runNumber,file)
	command = "bsub -q cmsexpress -o output_{1}.txt -e error_{1}.txt -J MinBias_2015_{1} automationMinBias.sh {0} {1}".format(runNumber,file)
        temp = os.popen(command)
	output = temp.read()
	jobid = output.split("<")[1].split(">")[0]
	command_wait = command_wait + " ended({0})".format(jobid)
	
os.system(command_wait+"\" automationMinBias_pede.sh {0} {1}".format(runNumber,tempNum))
#os.system('bash FindBanries.sh')
#os.system('bash automationMinBias_pede.sh')

#if os.path.exists("milleBinaryISN.dat"):
#    command = "rm milleBinaryISN.dat"
#    os.system(command)

os.chdir("..")

if sendEmailNotification:
	command = "echo \"New Alignment Updated for Run " + runNumber + "\" | mail -s \"New Prompt Alignment Update\" jcastle@cern.ch"
        os.system(command)

#if not upLoaded:
#    command = "tar -zcf " + directory + ".tgz " + directory
#    os.system(command)



#command = "rm -rf "+ directory + "/"
#os.system(command)

print 'Job Finished ' 
