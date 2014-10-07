#!/bin/bash
if [ $# -ne 6 ]
then 
  echo "Usage: ./makeDiJetIniSkim.sh <inputList> <sampleType> <outputFile> <outDir> <#> <justJtBool>"
  exit 1
fi

if [ $6 -eq 1 ]
then
  if [ $2 -eq 0 ]
  then
    tar -xzvf corrFilePbPb_20140429.tar.gz
  elif [ $2 -eq 1 ]
  then
    tar -xzvf corrFilePbPb_20140429.tar.gz
  elif [ $2 -eq 2 ]
  then
    tar -xzvf corrFilePP_20140430.tar.gz
  elif [ $2 -eq 3 ]
  then
    tar -xzvf corrFilePP_20140430.tar.gz
  fi
fi

echo | awk -v inputList=$1 -v sampleType=$2 -v outputFile=$3 -v num=$5 -v justJtBool=$6 '{print "./makeDiJetIniSkim.exe \""inputList"\" \""sampleType"\" \""outputFile"\" \""num"\" \""justJtBool"\""}' | bash
mv $3_$5*.root $4
rm *.root

echo "job done successfully"