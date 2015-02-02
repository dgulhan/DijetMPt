#!/bin/bash
if [ $# -ne 6 ]
then 
  echo "Usage: ./pmakeDiJetIniSkim.sh <inputList> <sampleType> <outName> <outDir> <justJtBool> <oldSample>"
  exit 1
fi

now="skimTreeJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $4
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetIniSkim.sh $now
cp ptCorrDirPbPb/*20140429*.tar.gz $now
# cp ptCorrDirPbPb/*20141220*.tar.gz $now
cp ptCorrDirPbPb/*20150118*.tar.gz $now
cp ptCorrDirPP/*20141105.tar.gz $now
cp ptCorrDirPP/*20141220.tar.gz $now
cp $1 $now

NAME="makeDiJetIniSkim.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe" `/net/hisrv0001/home/cfmcginn/FastJet/CMSSW_5_3_20/src/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

cp makeDiJetIniSkim.exe $now

cat pmakeDiJetIniSkim.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g" | sed "s@arg5@$5@g" | sed "s@arg6@$6@g" | sed "s@njobs@$len@g" > $now/pmakeDiJetIniSkim.condor
echo -=-
cat $now/pmakeDiJetIniSkim.condor
echo -=-
