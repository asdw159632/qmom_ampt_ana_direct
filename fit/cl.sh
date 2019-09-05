#!/bin/bash
pair=pi+pi+
cp ./pairs/$pair.h ./src/pair.h
make clean
make
if [ $? -ne "0" ]
then
	exit 0
fi
./analysis
