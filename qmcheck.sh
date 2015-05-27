#!/bin/bash
touch tempcheck1
touch tempcheck2
touch tempcheck3
touch tempcheck4
size1=0
size2=0
size3=0
size4=0
until [ "$size1" !=  "0"  -a "$size2" != "0" -a "$size3" != "0" -a "$size4" != "0" ] 
do
sleep 5
grep "Restart data updated" ./QMwork1/qm.log > tempcheck1
grep "Restart data updated" ./QMwork2/qm.log > tempcheck2
grep "Restart data updated" ./QMwork3/qm.log > tempcheck3
grep "Restart data updated" ./QMwork4/qm.log > tempcheck4
size1=`ls -l tempcheck1 | awk '{print $5}'`
size2=`ls -l tempcheck2 | awk '{print $5}'`
size3=`ls -l tempcheck3 | awk '{print $5}'`
size4=`ls -l tempcheck4 | awk '{print $5}'`
done
rm tempcheck*
echo "Two QM job finished"
