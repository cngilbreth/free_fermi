#!/bin/bash

PROG=../src/free_fermi_proj

betas="0.1 1.0 `seq 2.0 2.0 10.0`"
N=10

echo "#beta       E"
for beta in $betas; do
    echo -n "$beta   "
    $PROG E $N $beta
done

#echo ""
#echo "#beta       C"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG C $N $beta
#done

echo ""
echo "#beta       F"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F $N $beta
done

echo ""
for beta in $betas; do
    echo "#beta=$beta"
    $PROG nk $N $beta 20
    echo ""
done

#echo ""
#echo "#beta     E_spin"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG E_spin $N $beta
#done

#echo ""
#echo "#beta     F_spin"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG F_spin $N $beta
#done

