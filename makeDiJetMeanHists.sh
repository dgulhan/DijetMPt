#!/bin/bash

if [ $# -ne 5 ]
then 
  echo "Usage: ./makeDiJetMeanHists.sh <inputList> <sType> <outDir> <isHITrk> <#>"
  exit 1
fi

echo | awk -v inputList=$1 -v sType=$2 -v isHITrk=$4 -v num=$5 '{print "./makeDiJetMeanHists.exe \""inputList"\" \""sType"\" \""isHITrk"\" \""num"\""}' | bash

mv *MeanHists*.root $3
rm *.root 

echo "job done successfully"