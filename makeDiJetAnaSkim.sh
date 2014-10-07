#!/bin/bash

if [ $# -ne 6 ]
then 
  echo "Usage: ./makeDiJetAnaSkim.sh <inputList> <sType> <outDir> <#> <justJtBool> <isHITrk>"
  exit 1
fi

if [ $2 -eq 0 ] || [ $2 -eq 1 ]
then
    tar -xzvf corrFilePbPb_20140429.tar.gz
elif [ $2 -eq 2 ] || [ $2 -eq 3 ]
then
  if [ $7 -eq 0 ]
  then
      tar -xzvf corrFilePP_20140430.tar.gz
  elif [ $7 -eq 1 ]
  then
      tar -xzvf corrFilePP_HITrk_20140821.tar.gz
  fi
fi

echo | awk -v inputList=$1 -v sType=$2 -v num=$4 -v justJt=$5 -v isHITrk=$6 '{print "./makeDiJetAnaSkim.exe \""inputList"\" \""sType"\" \""num"\" \""justJt"\" \""isHITrk"\""}' | bash

rm cent*.root
rm hist*.root
rm eff*.root
rm fake*.root
mv *.root $3
rm *.root 

echo "job done successfully"