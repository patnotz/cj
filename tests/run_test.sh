#!/bin/bash

ESC_SEQ="\x1b["
COL_RESET=$ESC_SEQ"39;49;00m"
COL_RED=$ESC_SEQ"31;01m"
COL_GREEN=$ESC_SEQ"32;01m"
COL_YELLOW=$ESC_SEQ"33;01m"
COL_BLUE=$ESC_SEQ"34;01m"
COL_MAGENTA=$ESC_SEQ"35;01m"
COL_CYAN=$ESC_SEQ"36;01m"

echo -e $COL_CYAN"RUNNING REGRESSION TESTS"$COL_RESET
echo -e $COL_BLUE "$(date)" $COL_RESET
for file in *
do
  if [ "${file}" == "run_test.sh" -o  "${file}" == "run_test.sh~" -o "${file}" == "exodiff" ]
  then
    continue
  fi
  cd ./$file/
  ../../build/src/Cj
  ../exodiff -stat mesh.e ./gold/mesh.e  > test.log
  export ex_status="$?"
  if [ $ex_status == "0" ]
  then    
    echo -e $COL_GREEN "[ PASSED ] "$COL_RESET "${file}"
  elif [ $ex_status == "2" ]
  then
    echo -e $COL_YELLOW "[  DIFF  ] "$COL_RESET "${file}"    
  else
    echo -e $COL_RED "[ FAILED ] "$COL_RESET  "${file}"
  fi
  cd ../
done
