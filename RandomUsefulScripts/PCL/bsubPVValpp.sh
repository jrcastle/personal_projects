RUNNUM="$1"
eval `scram runtime -sh`
bsub -o log.log -q cmscaf1nd PVValidationSubmitter.csh /afs/cern.ch/work/j/jcastle/CMSSW_7_5_4/src/Alignment/OfflineValidation/test/InputSource_pp3.dat pp3.8T_PCL_Run$RUNNUM
bsub -o log2.log -q cmscaf1nd PVValidationSubmitter.csh /afs/cern.ch/work/j/jcastle/CMSSW_7_5_4/src/Alignment/OfflineValidation/test/InputSource_pp4.dat pp3.8T_ExpressGT_Run$RUNNUM