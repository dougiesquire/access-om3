#!/bin/bash

find MOM6 -type f -name "*.new" | while read -r file
do
    new=$file
    orig=${file%.new}
    fname=$(basename $orig)
    echo "Making patch for ${fname}"
    git diff --no-index --patch --output=patches/${fname}.patch ${orig} ${new}
done
