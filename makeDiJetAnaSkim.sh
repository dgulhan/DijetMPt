#!/bin/bash

CMSSWENV=/net/hidsk0001/d00/scratch/dgulhan/CMSSW_5_3_20/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc5_amd64_gcc462
cd $CMSSWENV
eval `scramv1 runtime -sh` 
cd -

if [ $# -ne 6 ]
then 
  echo "Usage: ./makeDiJetAnaSkim.sh <inputList> <sType> <outDir> <#> <justJtBool> <isHITrk>"
  exit 1
fi

if [ $2 -eq 0 ] || [ $2 -eq 1 ]
then
    tar -xzvf corrFilePbPb_20141219.tar.gz
   # tar -xzvf corrFilePbPb_20140429.tar.gz
   tar -xzvf swapCorrPbPbv5.tar.gz
#    tar -xzvf corrFilePbPb_20141209_HI.tar.gz
elif [ $2 -eq 2 ] || [ $2 -eq 3 ]
then
  if [ $6 -eq 0 ]
  then
      tar -xzvf corrFilePP_20141105.tar.gz
      tar -xzvf swapCorrPPv5.tar.gz
  elif [ $6 -eq 1 ]
  then
      tar -xzvf corrFilePP_HITrk_20141105.tar.gz
  fi
fi

echo | awk -v inputList=$1 -v sType=$2 -v num=$4 -v justJt=$5 -v isHITrk=$6 '{print "./makeDiJetAnaSkim.exe \""inputList"\" \""sType"\" \""num"\" \""justJt"\" \""isHITrk"\""}' | bash

rm cent*.root
rm hist*.root
rm eta*.root
rm eff*.root
rm fake*.root
rm secondary*.root
mv *Skim*.root $3
rm *.root

echo "job done successfully"