#!/bin/bash
#First argument corresponds to executing script
#Second argument corresponds to file listing files to execute script on

python /nfs/users/nfs_l/lbc/Codes/MyCodes/JobArrays/JobArrays.py ${1} ${2} ${LSB_JOBINDEX}
echo "python ${1} on line ${LSB_JOBINDEX} of ${2}"
