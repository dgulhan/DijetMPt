#!/bin/bash

if [ $# -ne 5 ]
then 
  echo "Usage: ./pmakeDiJetAnaSkim.sh <inputList> <sType> <outDir> <justJtBool> <isHITrk>"
  exit 1
fi

now="skimTreeJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $3
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetAnaSkim.sh $now
cp ptCorrDirPbPb/*.tar.gz $now
cp ptCorrDirPP/*.tar.gz $now
cp centWeightDir/centHist_eventSet_*.root $now
cp histWeightFile.root $now
cp $1 $now

NAME="makeDiJetAnaSkim.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe" `/net/hisrv0001/home/cfmcginn/FastJet/CMSSW_5_3_12_patch3/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
cp makeDiJetAnaSkim.exe $now

cat pmakeDiJetAnaSkim.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g" | sed "s@arg5@$5@g" | sed "s@njobs@$len@g" > $now/pmakeDiJetAnaSkim.condor
echo -=-
cat $now/pmakeDiJetAnaSkim.condor
echo -=-
