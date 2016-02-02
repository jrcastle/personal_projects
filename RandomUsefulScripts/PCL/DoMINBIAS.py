#!/usr/bin/env python
import sys
import subprocess
import os
from ROOT import TH1F, TFile, TH1
from subprocess import *
import time
import shlex
import stat

RUNNUM = sys.argv[1]
runNumber = sys.argv[2]
sendEmailNotification = False
start = time.time()
print 'Start New Job ' + str(time.asctime())

FileNames=open('TkMinBias'+runNumber,'r')
ListFiles = FileNames.readlines()
FileNames.close()

print 'Run Number ' + runNumber

# Template alignment cfg
TemplateCfg=open('alignment_mille_minbias.py','r')
TemplateCfgLines=TemplateCfg.readlines()
TemplateCfg.close()

# New alignment cfg that will have the DAS query files input
NewCfg=open('alignment_mille.py','w')

# Go through the template configuration file (alignment_BASE.py) and add the file names for the run into a new file (alignment_New.py)
for cfgLine in TemplateCfgLines:
    if 'INSERTFILENAME' in cfgLine:
        newline=''
        for fileName in ListFiles:
            newline = newline+cfgLine.replace('INSERTFILENAME',fileName[:-1])
    elif 'INSERT_NUM' in cfgLine:
    	newline = cfgLine.replace('INSERT_NUM',str(runNumber))
    elif 'RUN_NUMBER' in cfgLine:
	newline = cfgLine.replace('RUN_NUMBER',str(RUNNUM))
    else:
        newline=cfgLine
    NewCfg.write(newline)


NewCfg.close()

#Make directory for specific run and move files into there
directory = "MinBias_2015_" + runNumber

if os.path.exists(directory):
    command = "rm -rf "+directory
    os.system(command)

command = "mkdir " + directory
os.system(command)
command = "mv alignment_mille.py " + directory + "/."
os.system(command)
command = "cp TkMinBias"+runNumber+" " + directory + "/."
os.system(command)
os.getcwd()
os.chdir(directory)
os.getcwd()

#Run alignment
os.system('cmsRun alignment_mille.py')
end = time.time()

#if os.path.exists("milleBinaryISN.dat"):
#    command = "rm milleBinaryISN.dat"
#    os.system(command)

os.chdir("..")

if sendEmailNotification:
               command = "echo \"New Alignment Mille File for Run " + runNumber + "\" | mail -s \"New Prompt Alignment Update\" jcastle@cern.ch"
               print command
               os.system(command)
