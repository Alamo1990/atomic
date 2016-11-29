#!/bin/bash

make clean
make

rm -rf output
mkdir output

./lock_test 10

RESULT[0]=$(diff -s ../solution/0.txt output/0.txt)
RESULT[1]=$(diff -s ../solution/1.txt output/1.txt)
RESULT[2]=$(diff -s ../solution/2.txt output/2.txt)
RESULT[3]=$(diff -s ../solution/3.txt output/3.txt)
RESULT[4]=$(diff -s ../solution/4.txt output/4.txt)

RESULT[5]=$(diff -s ../solution/5.txt output/5.txt)
RESULT[6]=$(diff -s ../solution/6.txt output/6.txt)
RESULT[7]=$(diff -s ../solution/7.txt output/7.txt)
RESULT[8]=$(diff -s ../solution/8.txt output/8.txt)
RESULT[9]=$(diff -s ../solution/9.txt output/9.txt)

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


echo "-----------------------------LOCK-BASIC---------------------------"


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
