#This script uses all the output files from the picard CollectInsertSizeMetrics tool, assuming they have *_insert_size_txt extension (this can be changed below), in the directory.
##It will rely on the ShortInsert_Versus_LongInsert_v1.R Rscript, so please make sure you adapt the correpondent paths
# Lara Bossini-Castillo 08th Nov 2017 lbc@sanger.ac.uk
#!/bin/bash

files="*_insert_size_txt"

for i in ${files}
do
  sed '1,10d' ${i} > ${i}_2
  wait
  sed '/^$/d' ${i}_2 > ${i}_3 
  wait
  mv  ${i}_3 ${i}_2
done

wait

/software/R-3.2.1/bin/Rscript /nfs/users/nfs_l/lbc/Codes/MyCodes/ATAC-seq/ShortInsert_Versus_LongInsert_v1.R
