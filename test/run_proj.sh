#!/bin/bash

PROG=../src/free_fermi_proj

betas=`seq 1.0 1.0 10.0`
echo $betas

echo "#beta       E"
for beta in $betas; do
    echo -n "$beta   "
    $PROG E 10 $beta
done

#echo ""
#echo "#beta       C"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG C 10 $beta
#done

echo ""
echo "#beta       F"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F 10 $beta
done

#echo ""
#echo "#beta     E_spin"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG E_spin 10 $beta
#done

#echo ""
#echo "#beta     F_spin"
#for beta in $betas; do
#    echo -n "$beta   "
#    $PROG F_spin 10 $beta
#done

