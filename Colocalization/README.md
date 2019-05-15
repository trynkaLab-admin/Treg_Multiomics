This folder contains the scripts for colocalization analysis between Treg QTLs and immune variants. Since QTL analysis was done using GRCh38 and most of the summary statistics at the time of the analysis were in GRCh37, a large portion of the script relies on making the two compatible. The script was adapted from Kaur Alasoo colocWrapper (https://github.com/kauralasoo/colocWrapper). 

It requires a file listing all the GWAS studies that we want to test. The file should contain four columns, GWAS acronym, GWAS file name, type of GWAS and control_cases_ratio, such as:
CD	Crohns_disease_deLange_2017 Autoimmune	0.3028361
IBD	Inflammatory_bowel_disease_deLange_2017 Autoimmune	0.417666
UC	Ulcerative_cholitis_deLange_2017  Autoimmune	0.2689723

It requires two files listing all variants tested in the QTL study along with their MAF with GRCh37 and GRCh38 coordinates. Importantly, the SNP IDs should match, such as:
GRCh37.txt
1	751343	1_751343_T_A	T	A	SNP	0.138514
1	751756	1_751756_T_C	T	C	SNP	0.138514
1	752566	1_752566_G_A	G	A	SNP	0.152027

GRCh38.txt
1	815963	1_751343_T_A	T	A	SNP	0.138514
1	816376	1_751756_T_C	T	C	SNP	0.138514
1	817186	1_752566_G_A	G	A	SNP	0.152027
