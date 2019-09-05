#!/bin/csh

date

rm -fr ./tmp
mkdir ./tmp
rm -f result-all.root

@ i = 0
foreach list(`find ./result_root/ -name "*.root"`)
if ($i == 0) then
hadd ./tmp/tmp-all.root $list
mv ./tmp/tmp-all.root ./tmp/result-all.root
else
hadd ./tmp/tmp-all.root ./tmp/result-all.root $list
mv ./tmp/tmp-all.root ./tmp/result-all.root
endif
@ i += 1
end
mv ./tmp/result-all.root ./result_root/
rm -fr ~/tmp

echo "add successfully!"

date

#echo "input Number of Files"
#@ Nfiles = $<

#@ i = 0

#while($i<$Nfiles)

#  if ($i == 0) then
#  hadd tmp-all.root result-$i.root
#  mv tmp-all.root result-all.root
#  else
#  hadd tmp-all.root result-all.root result-$i.root
#  mv tmp-all.root result-all.root
#  endif

#  @ i += 1

#end
