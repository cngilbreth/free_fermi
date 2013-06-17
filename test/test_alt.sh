#!/bin/bash

# Test free_fermi_recurse

echo -n "Testing free_fermi_recurse ... "
./run_alt.sh > run_alt.txt
diff run_alt.txt run_alt1000.ref
if [ $? -ne 0 ]; then
    echo "Error: test failed!"
    echo "got                                                           expected"
    echo "---                                                           --------"
    diff -y run_alt.txt run_alt1000.ref
    exit 1
fi

echo "passed"
