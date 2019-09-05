#!/bin/bash
source ~/.bashrc
main=$(pwd ./)
#rawdatapath=/storage/staff3652/staff/szhang/ampt-deuteron-BES-DATA/AT197ZT79-AP197ZP79-200GeV-b0-4.416fm-ISOFT4-NTMAX150/2017050814/
pair=pi+pi+
rm -r AuAu
for i in AuAu
do
	rawdatapath=/media/zl/下载/HBT/ampt-initial-fluctuation-DATA/data-$i/
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
	make clean
	make -j4
	if [ $? -ne 0 ]
	then
		echo ""
		echo "Make Failed!"
		echo ""
		exit
	fi
	cd bin/
	./analysis $rawdatapath $pair
	cd $main
	mkdir $i
	cp -r result_* $i
done
#cd dpt-dy
#./cl.sh
#cd ..
cp -r AuAu/* ../myoutput/Triangle/
cp -r AuAu/* ../myoutput/Chain/
cp -r AuAu/* ../myoutput/nofluc/
cd ../myoutput/fit/
./cl.sh
exit
