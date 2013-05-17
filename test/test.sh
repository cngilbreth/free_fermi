#!/bin/bash

# Test free_fermi

echo -n "Testing free_fermi ... "
./run.sh > run.txt
diff run.txt run.ref
if [ $? -ne 0 ]; then
    echo "Error: test failed!"
    echo "got                                                           expected"
    echo "---                                                           --------"
    diff -y run.ref run.txt
    exit 1
fi

echo "passed"