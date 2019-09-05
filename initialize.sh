#!/bin/sh
source ~/.bashrc
source ~/.conenv.sh
rawdatapath=/storage/staff3652/staff/szhang/ampt-deuteron-BES-DATA/AT197ZT79-AP197ZP79-200GeV-b0-4.416fm-ISOFT4-NTMAX150/2017050814/
pair=pi+pi+
rm result_dat/*.dat
rm result_gif/*.gif
rm result_root/*.root
cp ./pairs/$pair.h ./src/pair.h
if [ $? -ne 0 ]
then
	echo ""
	echo "Copy pair.h Failed!"
	echo ""
	exit
fi
make
if [ $? -ne 0 ]
then
	echo ""
	echo "Make Failed!"
	echo ""
	exit
fi
cd bin/
#./analysis $rawdatalist $inputlist 1
./analysis $rawdatapath $pair
exit
