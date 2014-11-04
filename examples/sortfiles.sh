dir=$1

mkdir $dir/sortfiles

for FILE in `ls -v $dir`
do

#for FILE in $dir/*
#do
  read -a file <<< $(echo $FILE | tr "/" "\n")

  filee=${file[-1]}
  echo $filee
  sort $dir/$filee > $dir/sortfiles/$filee.txt

done

