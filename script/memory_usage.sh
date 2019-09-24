#!/bin/bash

if [ $# -eq 0 ]
then
    N=10
else
    N=$1
fi

LOCATION="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $LOCATION/../build

mkdir -p $LOCATION/../output

PROG="mesher_roi -p ../data/a.poly -u ../output/a -a "
DATE=$(date '+%Y-%m-%d_%H:%M:%S')
DIR="memory_usage_mesher_roi_$N_"$DATE
FILE="$DIR/memory_usage.txt"

mkdir -p $DIR

for i in $(seq 1 $N)
do
    valgrind --tool=massif --stacks=yes --massif-out-file=$DIR/massif.out.%p.$i ./$PROG $i
done

echo ""
touch $FILE
echo "Command : ${PROG}N" >> $FILE
echo "----------------------" >> $FILE
echo "N  | Peak memory usage" >> $FILE
echo "----------------------" >> $FILE
for i in $(seq 1 $N)
do
    RES=$(awk -f ../script/massif_analyser.awk $DIR/massif.out.*.$i)
    if [ $i -ge 10 ]
    then
        echo "$i | $RES" >> $FILE
    else
        echo "$i  | $RES" >> $FILE
    fi
done

cat $FILE