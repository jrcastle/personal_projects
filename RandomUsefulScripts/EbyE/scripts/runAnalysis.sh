DIR=$PWD
cd /home/j550c590/CMSSW_7_5_8_patch2/src
eval `scram runtime -sh`
cd $DIR

if [ -a V2Det.root ]
    then
    echo "Removing V2Det.root ..."
    rm V2Det.root
fi

if [ -a CastleEbyE.root ]
    then
    echo "Removing CastleEbyE.root ..."
    rm CastleEbyE.root
fi

if [ -a nohup.out ]
    then
    echo "Removing nohup.out"
    rm nohup.out
fi

root -l -b <<EOF
.x makeV2Det.C++
.q
EOF
root -l -b <<EOF
.x ReadTree_normDet.C++
.q
EOF