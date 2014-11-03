#!/bin/bash

if [ $# -ne 4 ]
then 
  echo "Usage: ./pmakeDiJetMeanHists.sh <inputList> <sType> <outDir> <isHITrk>"
  exit 1
fi

now="meanHistsJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $3
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetMeanHists.sh $now
cp $1 $now

NAME="makeDiJetMeanHists.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe" `/net/hisrv0001/home/cfmcginn/FastJet/CMSSW_5_3_12_patch3/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
cp makeDiJetMeanHists.exe $now

cat pmakeDiJetMeanHists.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g" | sed "s@njobs@$len@g" > $now/pmakeDiJetMeanHists.condor
echo -=-
cat $now/pmakeDiJetMeanHists.condor
echo -=-
