#!/bin/bash

platform=`uname`

if [[ $platform == 'Linux' ]]; then
    base64cmd="base64 -w0"
elif [[ $platform == 'Darwin' ]]; then
    base64cmd=base64
else
    echo "Unsupported OS - must be Mac or Linux"
    exit 1
fi
    
fname=$1

for img in `egrep -o "[A-Za-z0-9\-]+\.png" $fname`; do
    img=${img%.*}
    echo -n "data:image/png;base64," > $img.base64
    $base64cmd < $img.png >> $img.base64
    sed -e "s/$img\.png/$(sed 's:/:\\/:g' $img.base64)/" $fname > $fname.embed
    mv $fname.embed $fname
done



