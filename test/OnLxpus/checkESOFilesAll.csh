#!/bin/tcsh
if ( $2 == "" ) then
	echo ">> [INFO] Please input the data forder list."
	echo ">>        Ex: ./checkESOFilesAll.csh [data_card] [server]"
        echo ">>        [data_card]:"
	echo '>>        ................................................................'	
	echo '>>        :    In data_card you should input like as following           :'
	echo '>>        :    Ex: FUll_Path;Jobs_numbers                                :'
	echo '>>        :        /dpm/grid.sinica.edu.tw/.../...;1004                  :'
	echo '>>        ................................................................'
        echo ">>        [server]: cmslpc or lxplus"
	exit	
endif
if ( ! ( -e $1 ) ) then
	echo ">> [WARING] $1 not found, please check."
	exit
endif

if ( $2 != 'cmslpc' && $2 != 'lxplus' ) then
	echo ">> [ERROR] $2 not found server options: cmslpc or lxplus."
        exit
endif

rm -f check_log/status.txt 
touch check_log/status.txt 
set datasets=`cat "$1"`
foreach data($datasets)
	set name=`echo $data | awk -F ";" '{print $1}'`
	set size=`echo $data | awk -F ";" '{print $2}'`
        source checkESOFiles.csh $name $size $2 | tee -a check_log/status.txt
end
set finish=`cat check_log/status.txt | grep '100.0%'    | wc -l`
set done=`  cat check_log/status.txt | grep 'DONE'      | wc -l`
set dup=`   cat check_log/status.txt | grep 'Duplicate' | wc -l`
set all=`cat $1 | wc -l `
set effall=`echo "scale=1; $finish*100/$all" | bc`
echo "" | tee -a check_log/status.txt
echo ">@ ============================================================== @<" | tee -a check_log/status.txt
echo ">> [INFO] All Jobs: $all"                                             | tee -a check_log/status.txt
echo ">> [INFO] Finished: $finish (Done:$done / Duplicate:$dup)"            | tee -a check_log/status.txt
echo '>> [INFO] Total Status: '$effall'%'                                   | tee -a check_log/status.txt
echo ">@ ============================================================== @<" | tee -a check_log/status.txt
echo ">> [INFO] The error massage will store in the check_log"

