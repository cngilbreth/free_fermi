#!/bin/bash

PROG=../src/free_fermi

betas=`seq 1.0 1.0 10.0`
echo $betas

echo "#beta       H"
for beta in $betas; do
    echo -n "$beta   "
    $PROG H 10 10 $beta
done

echo ""
echo "#beta       C"
for beta in $betas; do
    echo -n "$beta   "
    $PROG C 10 10 $beta
done

echo ""
echo "#beta       F"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F 10 10 $beta
done

echo ""
echo "#beta     H_spin"
for beta in $betas; do
    echo -n "$beta   "
    $PROG H_spin 10 10 $beta
done

echo ""
echo "#beta     F_spin"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F_spin 10 10 $beta
done

