#!/bin/bash

if [ $# -ne 5 ]
then 
  echo "Usage: ./pmakeDiJetAnaSkim_MIXING.sh <inputList> <sType> <outName> <outDir> <justJtBool>"
  exit 1
fi

now="skimTreeJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $4
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetAnaSkim_MIXING.sh $now
cp $1 $now

NAME="makeDiJetAnaSkim_MIXING.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe" `/net/hisrv0001/home/cfmcginn/FastJet/CMSSW_5_3_12_patch3/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
cp makeDiJetAnaSkim_MIXING.exe $now

cat pmakeDiJetAnaSkim_MIXING.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g" | sed "s@arg5@$5@g" | sed "s@njobs@$len@g" > $now/pmakeDiJetAnaSkim_MIXING.condor
echo -=-
cat $now/pmakeDiJetAnaSkim_MIXING.condor
echo -=-
