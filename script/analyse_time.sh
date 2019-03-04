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