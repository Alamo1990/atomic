#!/bin/bash

cd build
make clean
make build

rm -rf build/output
mkdir build/output

./build/lock_test 10

RESULT[0]=$(diff -s solution/0.txt build/output/0.txt)
RESULT[1]=$(diff -s solution/1.txt build/output/1.txt)
RESULT[2]=$(diff -s solution/2.txt build/output/2.txt)
RESULT[3]=$(diff -s solution/3.txt build/output/3.txt)
RESULT[4]=$(diff -s solution/4.txt build/output/4.txt)

RESULT[5]=$(diff -s solution/5.txt build/output/5.txt)
RESULT[6]=$(diff -s solution/6.txt build/output/6.txt)
RESULT[7]=$(diff -s solution/7.txt build/output/7.txt)
RESULT[8]=$(diff -s solution/8.txt build/output/8.txt)
RESULT[9]=$(diff -s solution/9.txt build/output/9.txt)

for i in 0 1 2 3 4 5 6 7 8 9
do
    AUX=${RESULT[i]}
    if [ "${AUX:(-9)}" == "identical" ] || [ "${AUX:(-9)}" == "idénticos" ]
    then
        RESULT[i]="Correct!!! =-D"
    else
        RESULT[i]="fail :("
    fi

done


echo "-----------------------------LOCK---------------------------"


echo "Image 0:  ${RESULT[0]}"
echo
echo "Image 1:  ${RESULT[1]}"
echo
echo "Image 2:  ${RESULT[2]}"
echo
echo "Image 3:  ${RESULT[3]}"
echo
echo "Image 4:  ${RESULT[4]}"
echo
echo "Image 5:  ${RESULT[5]}"
echo
echo "Image 6:  ${RESULT[6]}"
echo
echo "Image 7:  ${RESULT[7]}"
echo
echo "Image 8:  ${RESULT[8]}"
echo
echo "Image 9:  ${RESULT[9]}"
echo
echo
