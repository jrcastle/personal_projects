RUNNUM="$1"
rm -r MinBias_2015
scp -r jcastle@lxplus.cern.ch:/afs/cern.ch/user/j/jcastle/public/HeavyIonRun2015_PseudoPCL/pp3.8T_PCL_Alignment/Results$RUNNUM/MinBias_2015 .
echo "DONE!"