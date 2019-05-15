#### This is a Genotype refinement and Imputation + Phasing BEAGLE 4.1 pipeline that has been used for BLUEPRINT Phase 2 data.

### Pipeline details...

**Prepared by:** Kousik Kundu \
**Date:** August, 2017 \
**Purpose:** Genotype refinement and Imputation+phasing \
**Which project this pipeline has been applied to:** BLUEPRINT Phase 2 \
**Software:** BEAGLE 4.1 (```/nfs/team151/software/beagle_4/beagle.21Jan17.6cc.jar```)

--------------------------------------------------------------------
## Genotype refinement and Imputation + Phasing

***Important note:*** When using gl= input data for genotype refinement, the output is unphased genotypes and you will need to run Beagle a second time (with the gt= parameter) to phase the resulting genotypes.  The second run (with gt=) should take much less time than the first run. Beagle estimates the phase in the first run (with the gl = parameter) but reports the genotypes as unphased genotypes because you will get more accurate phase estimates with the two-stage strategy, in which genotypes are estimated in the first run and phased in the second run (with the gt = parameter).


### Beagle Pipeline 

I used Beagle pipeline (created by Petr Danecek) to run genotype refinement and Imputation + Phasing.

Github source for vt-runner: \
https://github.com/VertebrateResequencing/vr-runner/

In more detail, please read in the GitHub help page \
**vt-runner git vertion** - commit a43a14ad3d80b987785d61a7d603726807d9e492

```
# Get the code
git clone git://github.com/VertebrateResequencing/vr-runner.git

# Create a config file, edit, and run
cd /nfs/team151/software/vr-runner/
run-st-mpileup +sampleconf > my.conf
```

--------------------------------------------------------------------

### Steps for genotype refinement, imputation + phasing

1. [Setting environment variables](#STEP1)
2. [Creating and editing myGL.conf file](#STEP2)
3. [Input data and reference panel](#STEP3)
4. [RUN the pipeline for genotype refinement](#STEP4)
5. [Indel Alignement](#STEP5) 
6. [Editing myGT.conf file](#STEP6) 
7. [Input data and reference panel](#STEP7) 
8. [RUN the pipeline again for Imputation](#STEP8) 
9. [Second round phasing](#STEP9)
10. [Genotype concordance check](#STEP10) 

***Please note:*** [Step 5](#STEP5) should be done before running Beagle. We realized it after running genotype refinement step, and since genotype refinement step is very lengthy process, we decided to run this step (indel aligenment) after genotype refinement and before inputation+phasing. 

--------------------------------------------------------------------

<a name="STEP1"></a>
***STEP 1: Setting environment variables***

Set environmental variable in this file (if not present, create): ~/.vrw/runners
```
cat ~/.vrw/runners

#!/bin/bash
export PATH=/nfs/team151/software/vr-runner/scripts:$PATH
export PERL5LIB=/nfs/team151/software/vr-runner/modules:$PERL5LIB
```

The vr-wrapper script can be used to set environment variables without having to change user's profile
I also add this variables in my .bashrc file

--------------------------------------------------------------------
<a name="STEP2"></a>
***STEP 2: Creating and editing myGL.conf file***

```
run-beagle +sampleconf > myGL.conf

EDIT: IMPORTANT
Used the myGL.conf for Genotype refinement. If you use again (example for imputation/phasing, explaind in step 6) check below mentioned line carefully -

• beagle_jar => 
• beagle_args =>
• known_vcf =>
• input_tag => 'gl' - should be gl
• region =>
• limits => normally use runtime => 24*60, which submit jobs in 'long' queue. However, this can be increased 72*60 for 'basement' queue.

```
--------------------------------------------------------------------
<a name="STEP3"></a>
***STEP 3: Input data and reference panel***

```
Input data: \
Before running the imputation, we only used the VSQLOD filter, SNPGap, and IndelGap filters. \
I run this all 200 BLUEPRINT sampels.

Input data path: \
/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/VCF_QC/GATK_VQSR/{CHROM}.vqsr.PASS.vcf.gz

Reference panel path: \
/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf/{CHROM}.bcf
```
--------------------------------------------------------------------
<a name="STEP4"></a>
***STEP 4: RUN the pipeline for genotype refinement***

Run it on SCREEN

 • Launch a new screen by executing the command screen (open terminal and type 'screen') \
 • Type enter and run whatever commands you like. (In our case it will be the run-st-mpileup command) \
 • Detach from the screen by pressing ctrl+a followed by d (Ctrl + a +d) \
 • At any time one can reattach to a running screen by executing 'screen -r' 

If the farm has multiple head nodes, note that the screen lives on a concrete computer and therefore one has to log into the exact same machine to access it. One can simultaneously run as many screens as necessary.

Note that it is not necessary (nor recommended) to submit the runner scripts to the farm as they take little CPU and memory. Nonetheless, if the pipeline is run this way, the exit status is 111 when successfully completed.

```
RUN:
export LSB_DEFAULTGROUP=hematopoiesis
vr-wrapper ~/.vrw/runners /nfs/team151/software/vr-runner/scripts/run-beagle +config myGL.conf +loop 900 +retries 10 +maxjobs 500 -i Chr_GL/{CHROM}.vqsr.PASS.vcf.gz  -o Beagle_output_GL 2>&1 | tee -a full-log.GL.txt
```

Runtime script doesn't work properly in the pipeline. The script is here if you want to use it, it's quite handy sometimes in general:
```/software/vertres/codebase/scripts//runtime-stats```

**Runtime:** 16 days (with modelscale parameter set to **2.0**. Used beagle version is ```/nfs/team151/software/beagle_4/beagle.21Jan17.6cc.jar```)


**Modelscale parameter**
*modelscale=[positive number] specifies the model scale parameter when sampling haplotypes for unrelated individuals (default: modelscale=0.8). Increasing the modelscale parameter will trade reduced phasing accuracy for reduced run-time. However, when estimating posterior probabilities from genotype likelihood data, increasing the modelscale parameter could improve both accuracy and run-time.*

The Beagle developer (Brian Browning <browning@u.washington.edu>) suggested any value 1.2 to 2.5 should be fine. However, after 2.0 the performance would be monotonically reduced. He also said - \
*"The model scale parameter controls the number of states per level for the DAG model by controlling the threshold used to merge DAG nodes (our earliest Beagle papers contain some discussion of this).  In general increasing the model scale parameter (i.e. decreasing the number of model states) decreases phasing accuracy. However I have had a report from user who found a small increase in the model scale parameter (probably to some value between 1.0 and 2.0) could increase accuracy of genotype refinement."*

Note that in default modelscale parameter, the will take ages (guess 6 months or more :P).

**Pipeline progress**

• First chunks all the chromosome in a smaller regions. That can be decided in buffer_nsites and chunk_nsites in my.conf (I used a 1000 + 3000).

• Then run beagle on the chunk files separately.

• Finally merges the chunk files for each chromosome, and produces 22 VCF files for 22 chromosomes.

• Important think to note that Beagle 'gl' parameter only consider the intersection variants in the input and reference panel data. So in this case BLUEPRINT specific data that are not present is the reference panel will be ignored.

• The pipeline adds those BLUEPRINT specific variants to the Beagle output, which means the Beagle output files contain BLUEPRINT specific variants that are not refined (no AR2 value).

• NO phasing has been done is this stage.



**While RUNNING** \
There might be several LSF uses take place, like jobs go in UNKWN, USUSP, ZOMBIE status. One possible way to the job that has problem; the pipeline find it automatically and resubmit the particular job. However sometimes, that might not happen, and pipeline may looses track of some jobs. The best think would be to wait till there is no jobs running and left in the queue. Then kill the Beagle_output/.jobs folder and re-run the pipeline again. The pipeline is smart enough to see what are finished and what are left.

--------------------------------------------------------------------
<a name="STEP5"></a>
***STEP 5: Indel Alignement***

Ideally, we should have done this before genotype refinement or may be even before GATK VQSR. This is an important step, we realize later stage where we observed we are getting much less number of INDELs in compare to  reference database. Moreover, there are lots of BLUEPRINT specific INDELs that are not present in reference panel. We investigated and found many example like below -

Different Ref/Alt of an INDEL (same chr position) are present in BP and UK10K+1kGPIII. Hence, in the imputation process, beagle ignored them, and eventually they were considered as BP specific, where the same chromosome position is present in reference panel with difference ref and alt.
```
Example:
In BP data:  2_233138541_TC_TCC
In Ref data: 2_233138541_T_TC
```

To address this issue, we have used VT (https://github.com/atks/vt)

```
Input data: Beagle output data from Step 4
Output data: Save in "INDEL_left_alignment" directory.

RUN:
/nfs/team151/software/vt/vt_v0.5/vt normalize \
input_file <Beagle_output_GL/20.vcf.gz>
-r /lustre/scratch114/resources/ref/Homo_sapiens/1000Genomes_hs37d5/hs37d5.fa \
-o output_file <INDEL_left_alignment/20.indel_aligned.vcf.gz>
```
--------------------------------------------------------------------

<a name="STEP6"></a>
### Imputation and phasing Process:

***STEP 6: Editing myGT.conf file***

```
Copy from my.conf file.
EDIT:
• Set model scale parameter as default
• input_tag => 'gt' - should be gt
• Rest all are same as myGL.conf

```
--------------------------------------------------------------------
<a name="STEP7"></a>
***STEP 7: Input data and reference panel***

```
Input data: Beagle output data after indel alignment from Step 5

Reference panel path: Same as Step 4 (/lustre/scratch116/vr/ref/human/GRCh37/imputation/uk10k+1000g-phase3/vcf/{CHROM}.bcf)
```
--------------------------------------------------------------------
<a name="STEP8"></a>
***STEP 8: RUN the pipeline again for Imputation***

Run it on SCREEN
open terminal and type 'screen'
detach from screen (Ctrl + a +d)
return to screen (Ctrl + r)

```
RUN:
export LSB_DEFAULTGROUP=hematopoiesis
vr-wrapper ~/.vrw/runners /nfs/team151/software/vr-runner/scripts/run-beagle +config myGT.conf +loop 900 +retries 10 +maxjobs 500 -i Chr_GT/{CHROM}.indel_aligned.vcf.gz   -o Beagle_output_GT_indel_aligned 2>&1 | tee -a ull-log.GT.indel_aligned.txt
```

**Runtime:** 2 days

**Pipeline progress**

• First chunks all the chromosome in a smaller regions. That can be decided in buffer_nsites and chunk_nsites in my.conf (I used a 1000 + 3000). Same as before.
	
• Then run beagle on the chunk files separately.

• Finally merges the chunk files for each chromosome, and produces 22 VCF files for 22 chromosomes.

• Important think to note that Beagle output with the 'gt' parameter does not contain BLUEPRINT specific data that are not present is the reference panel.

• The pipeline adds those BLUEPRINT specific variants to the Beagle output, which means the Beagle output files contain BLUEPRINT specific variants that are not PHASED. However, I have observed some of the SNPs are phased for some samples (wired; need to check), while input data is unphased. This is not Beagle problem. There might be something in the pipeline that needs to be checked. It doesn't affect us, as I filtered out all Blueprint specific variants for the downstream analyses.


--------------------------------------------------------------------
<a name="STEP9"></a>
***STEP 9: Second round phasing***

Because of the unphased BLUEPRINT specific SNPs are observed in final imputed files, 2nd time phasing can be run.
This time do not use "Ref=" parameter in Beagle; only use "gt"

**NOTE**: We haven't run this step, as we ignore all the BLUEPRINT specific variants that are not present is reference panel (UK10K+1KGPIII)

--------------------------------------------------------------------
<a name="STEP10"></a>
***STEP 10: Genotype concordance check***

We used GenotypeConcordance from GATK.

```
 java -jar GenomeAnalysisTK.jar \
   -T GenotypeConcordance \
   -R [hg19 reference].fasta \
   -eval [after genotype refinement].vcf \
   -comp [before genotype refinement].vcf \
   -o output_concordance.out
```
For more details, please see [GATK tool documentation](https://software.broadinstitute.org/gatk/gatkdocs/3.5-0/org_broadinstitute_gatk_tools_walkers_variantutils_GenotypeConcordance.php)

• Check the genotype concordance before and after genotype refinement. It gives us an idea that how much genotypes has been changed after genotype refinement. In our case around 1.5 to 2% genotypes were changed with NRD 6.3%.

• Check the genotype concordance before and after imputation and phasing. It gives us an insight, if there is any significant changes are observed, which is not expected. In our case we achieved very good concordance 99.9% (NRD 0.47%).

--------------------------------------------------------------------
