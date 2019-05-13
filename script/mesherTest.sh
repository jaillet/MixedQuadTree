#!/usr/bin/env bash

LOCATION="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $LOCATION/../build

PROG="mesher_roi"
DATE=$(date '+%Y-%m-%d_%H:%M:%S')
DIR="analyse_mesher_"$DATE
MAXRL=13


mkdir -p $DIR

for OPTION in s a
do
	echo "Launch mesher_roi on a.poly with option : -$OPTION $MAXRL"
	mkdir -p $DIR/${OPTION}_${MAXRL}

    for TRY in {1..10}
    do
    	echo "Try number $TRY"
        ./$PROG -p ../data/a.poly -$OPTION $MAXRL > $DIR/${OPTION}_${MAXRL}/test_number_${TRY}.txt
      
    done
done