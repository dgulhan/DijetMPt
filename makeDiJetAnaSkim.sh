if [ $# -ne 7 ]
then 
  echo "Usage: ./makeDiJetAnaSkim.sh <inputList> <sType> <outputFile> <outDir> <#> <justJtBool> <isHITrk>"
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

echo | awk -v inputList=$1 -v sType=$2 -v outputFile=$3 -v num=$5 -v justJt=$6 -v isHITrk=$7 '{print "./makeDiJetAnaSkim.exe \""inputList"\" \""sType"\" \""outputFile"\" \""num"\" \""justJt"\" \""isHITrk"\""}' | bash
mv $3_$5.root $4
rm *.root

echo "job done successfully"