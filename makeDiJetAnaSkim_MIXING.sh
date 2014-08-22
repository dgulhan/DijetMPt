if [ $# -ne 6 ]
then 
  echo "Usage: ./makeDiJetAnaSkim_MIXING.sh <inputList> <sType> <outputFile> <outDir> <#> <justJtBool>"
  exit 1
fi

echo | awk -v inputList=$1 -v sType=$2 -v outputFile=$3 -v num=$5 -v justJt=$6 '{print "./makeDiJetAnaSkim_MIXING.exe \""inputList"\" \""sType"\" \""outputFile"\" \""num"\" \""justJt"\""}' | bash
mv $3_$5.root $4
rm *.root

echo "job done successfully"