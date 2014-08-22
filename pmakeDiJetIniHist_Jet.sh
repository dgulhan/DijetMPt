if [ $# -ne 4 ]
then 
  echo "Usage: ./makeDiJetIniHist_Jet.sh <inputList> <outName> <outDir> <isPbPb>"
  exit 1
fi

now="iniHistJob_$(date +"%m_%d_%Y__%H_%M_%S")"
mkdir $now
mkdir -p $3
len=`wc -l $1 | awk '{print $1}'`
cp makeDiJetIniHist_Jet.sh $now
cp $1 $now

NAME="makeDiJetIniHist_Jet.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"
cp makeDiJetIniHist_Jet.exe $now

cat pmakeDiJetIniHist_Jet.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@arg1@$1@g" | sed "s@arg2@$2@g" | sed "s@arg3@$3@g" | sed "s@arg4@$4@g" | sed "s@njobs@$len@g" > $now/pmakeDiJetIniHist_Jet.condor
echo -=-
cat $now/pmakeDiJetIniHist_Jet.condor
echo -=-
