#!/usr/bin/env bash

LOCATION="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $LOCATION/../build

PROG="meanTime "
DATE=$(date '+%Y-%m-%d_%H:%M:%S')
DIR="analyse_time_"$DATE
FILE="$DIR/analyse_time.txt"

mkdir -p $DIR

for SIZE in 100 1000 10000 100000 1000000 10000000 100000000 1000000000
do
    for THREAD in 2 4 8 16
    do
        RES=$(./$PROG $SIZE $THREAD)

        echo "$RES" > time_size_${SIZE}_thread_${THREAD}.txt
    done
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