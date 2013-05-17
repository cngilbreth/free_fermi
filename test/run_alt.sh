#!/bin/bash

PROG=../altsrc/free_fermi_recurse

betas="0.1 1.0 2.0 4.0 6.0 8.0 10.0"
N=10

echo "#beta       E"
for beta in $betas; do
    echo -n "$beta   "
    $PROG E $N $beta
done

echo ""
echo "#beta       C"
for beta in $betas; do
    echo -n "$beta   "
    $PROG C $N $beta
done

echo ""
echo "#beta       F"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F $N $beta
done

echo ""
echo "#beta     E_spin"
for beta in $betas; do
    echo -n "$beta   "
    $PROG E_spin $N $beta
done

echo ""
echo "#beta     F_spin"
for beta in $betas; do
    echo -n "$beta   "
    $PROG F_spin $N $beta
done

