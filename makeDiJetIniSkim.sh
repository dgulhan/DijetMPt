#!/bin/bash
if [ $# -ne 7 ]
then 
  echo "Usage: ./makeDiJetIniSkim.sh <inputList> <sampleType> <outputFile> <outDir> <#> <justJtBool> <oldSample>"
  exit 1
fi

if [ $2 -eq 0 ]
then
    tar -xzvf corrFilePbPb_20140429.tar.gz
    tar -xzvf FFJEC_correction_PF_akVs_PbPb_20150118.tar.gz
    tar -xzvf residualcorr_akVs_PbPb_20150118.tar.gz
    tar -xzvf effCorrNPF_PbPb_20150118.tar.gz
elif [ $2 -eq 1 ]
then
    tar -xzvf corrFilePbPb_20140429.tar.gz
    tar -xzvf FFJEC_correction_PF_akVs_PbPb_20150118.tar.gz
    tar -xzvf residualcorr_akVs_PbPb_20150118.tar.gz
    tar -xzvf effCorrNPF_PbPb_20150118.tar.gz
elif [ $2 -eq 2 ]
then
    tar -xzvf corrFilePP_20140430.tar.gz
elif [ $2 -eq 3 ]
then
    tar -xzvf corrFilePP_20140430.tar.gz
fi


echo | awk -v inputList=$1 -v sampleType=$2 -v outputFile=$3 -v num=$5 -v justJtBool=$6 -v oldSample=$7 '{print "./makeDiJetIniSkim.exe \""inputList"\" \""sampleType"\" \""outputFile"\" \""num"\" \""justJtBool"\" \""oldSample"\""}' | bash
mv $3_$5*.root $4
rm *.root

echo "job done successfully"