#!/bin/sh
RUNNUM="$1"
NUMFILES="$2"
cd /afs/cern.ch/work/j/jcastle/CMSSW_7_5_5_patch1/src/
eval `scram runtime -sh`
cd /afs/cern.ch/work/j/jcastle/CMSSW_7_5_5_patch1/src/PbPb3.8T_2015_Alignment/Results$RUNNUM/

bash findBinaries.sh $RUNNUM $NUMFILES
python DoMINBIAS_pede.py $RUNNUM
mv DQM_plots.py /afs/cern.ch/work/j/jcastle/CMSSW_7_5_5_patch1/src/PbPb3.8T_2015_Alignment/Results$RUNNUM/MinBias_2015/
cd /afs/cern.ch/work/j/jcastle/CMSSW_7_5_5_patch1/src/PbPb3.8T_2015_Alignment/Results$RUNNUM/MinBias_2015/
python DQM_plots.py
