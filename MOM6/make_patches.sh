#!/bin/bash

sed -i '/add_patched_source/d' CMakeLists.txt
find MOM6 -type f -name "*.new" | while read -r file
do
    new=$file
    orig=${file%.new}
    fname=$(basename $orig)

    echo "Making patch for ${fname}"
    sed -i "\|$orig|d" CMakeLists.txt
    git diff --no-index --patch --output="patches/${fname}.patch" ${orig} ${new}
    echo "add_patched_source(OM3_mom6 ${orig})" >> CMakeLists.txt
done
