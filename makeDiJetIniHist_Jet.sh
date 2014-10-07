#!/bin/bash

if [ $# -ne 5 ]
then 
  echo "Usage: ./makeDiJetIniHist_Jet.sh <inputList> <outFileName> <outDir> <isPbPb> <isMonteCarlo>"
  exit 1
fi

echo | awk -v inputList=$1 -v outFileName=$2 -v isPbPb=$4 -v isMonteCarlo=$5 '{print "./makeDiJetIniHist_Jet.exe \""inputList"\" \""outFileName"\" \""isPbPb"\" \""isMonteCarlo"\""}' | bash
mv $2.root $3

echo "job done successfully"