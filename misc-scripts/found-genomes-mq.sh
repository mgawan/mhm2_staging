#!/bin/bash

r=$1
prefix=$2
echo "bin1 name refdepth genfrac"
for b in `seq 0 10000`; do
    f=mq-${prefix}_$b/summary/TXT/Genome_fraction.txt
    if [ ! -f "$f" ]; then
        continue
    fi
    x=`grep -v "-" $f`
    echo $x
    bin=`echo $x|cut -d' ' -f2|cut -d '_' -f2`
    echo -n "$bin"
    for i in `echo $x|cut -d' ' -f3,4`; do
        y=`grep $i $r|lastcol.py`
        echo -n " $i $y"
    done
    echo ""
done
