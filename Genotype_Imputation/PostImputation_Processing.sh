#!/usr/bin/env bash
#First argument is an interger supplying chromosome number
#Second argument is an interger [0-1] establishing MAF threshold

#Add tags to vcf
#bcftools +fill-tags $1.vcf.gz -- -t AN,AC,AF,HWE > $1_Test.vcf
#wait

#Filter by AR2 and HWE
#bcftools view ../$1_Test.vcf -M 2 -c 5 -i 'AR2>=0.8 && HWE>=1e-3' -Oz -o $1_AR2_HWE_filtered.vcf.gz
#wait

##Add GRCh37 position to SNP name
#gunzip $1_AR2_HWE_filtered.vcf.gz
#wait
#vim -c"31,\$s/\(\d\+\)\s\+\(\d\+\)\s\+\(\.\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(.\+\)/chr\1\t\2\tchr\1_\2_\4_\5\t\4\t\5\t\6/e|wq" $1_AR2_HWE_filtered.vcf
#wait
#vim -c"31,\$s/:\S\+\s\+/\t/g|wq" $1_AR2_HWE_filtered.vcf
#wait
#vim -c"31,\$s/:\S\+$//e|wq" $1_AR2_HWE_filtered.vcf
#wait

###Filtering by MAF
#plink --memory 4000 --vcf $1_AR2_HWE_filtered.vcf --maf $2 --make-bed --out $1_AR2_HWE_filtered_MAF_$2 
#wait

###Lifting over GRCh38
#awk -F"\t" '{print "chr" $1, $4, $4+1, $2}' $1_AR2_HWE_filtered_MAF_$2.bim > $1_AR2_HWE_filtered_MAF_$2.BED
#wait
#vim -c"%s/\s\+/\t/g|wq" $1_AR2_HWE_filtered_MAF_$2.BED 
#wait
#/software/pathogen/external/apps/usr/bin/liftOver $1_AR2_HWE_filtered_MAF_$2.BED /lustre/scratch117/cellgen/teamtrynka/lara/Temp_Resources/hg19ToHg38.over.chain.gz $1_AR2_HWE_filtered_MAF_$2_GRCh38.BED $1_AR2_HWE_filtered_MAF_$2_GRCh38.unmapped
#wait
#cut -f4 $1_AR2_HWE_filtered_MAF_$2_GRCh38.unmapped > $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove
#wait
#grep -v Deleted $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove > $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove2
#wait
#mv $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove2 $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove
#wait
#awk -F"\t" '{print $4, $2}' $1_AR2_HWE_filtered_MAF_$2_GRCh38.BED > $1_AR2_HWE_filtered_MAF_$2_GRCh38.UpdatePosition
#wait
#plink --memory 4000 --bfile $1_AR2_HWE_filtered_MAF_$2 --exclude $1_AR2_HWE_filtered_MAF_$2_GRCh38.2Remove --update-map $1_AR2_HWE_filtered_MAF_$2_GRCh38.UpdatePosition --make-bed --out $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test 
#wait
###Remove alleles from name and variants with the same position
#cut -f4 $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test.bim | uniq -d > $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test.bim2
#wait
#grep -Fwf $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test.bim2 $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test.bim > $1_SNPs2exclude.txt
#wait
#plink --memory 4000 --bfile $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test --exclude $1_SNPs2exclude.txt --make-bed --out $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos
#wait
#cp $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos.bim $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos.bim2
#wait
#vim -c"%s/\(\S\+\)\s\+\(\S\+_\d\+_\)\(\S\+_\S\+\)\s\+\(.\+\)/\1\t\2\t\4/e|wq" $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos.bim2
#wait
#mv $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos.bim2 $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos.bim
#wait
##Force reference allele
plink --memory 4000 --bfile $1_AR2_HWE_filtered_MAF_$2_GRCh38_Test_NoDupPos --exclude /lustre/scratch117/cellgen/teamtrynka/lara/Temp_Resources/BLUEPRINT_Beagle_Ref/UK10K_maf0.0001_chr$1_GRCh38.TriAllelic --a1-allele /lustre/scratch117/cellgen/teamtrynka/lara/Temp_Resources/BLUEPRINT_Beagle_Ref/UK10K_maf0.0001_chr$1_GRCh38.AltAlleles2 --recode vcf --out $1_AR2_HWE_filtered_MAF_$2_GRCh38_Setted
wait
##Add alleles to SNP names
vim -c"8,\$s/\(\d\+\)\s\+\(\d\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(\S\+\)\s\+\(.\+\)/chr\1\t\2\t\3\4_\5\t\4\t\5\t\6/e|wq"  $1_AR2_HWE_filtered_MAF_$2_GRCh38_Setted.vcf 
wait

echo "Setting of chromosome $1 Done!"
