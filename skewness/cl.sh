#!/bin/bash
make clean
make
if [ $? -ne "0" ]
then
	exit 0
fi
./analysis
