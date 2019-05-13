#!/usr/bin/env bash

LOCATION="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $LOCATION/../build

PROG="mesher_roi"
MAXRL=15

for OPTION in s a
do
	echo "Launch mesher_roi on a.poly with option : -$OPTION $MAXRL"
    ./$PROG -p ../data/a.poly -$OPTION $MAXRL
done