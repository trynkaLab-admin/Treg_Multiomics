#!/usr/bin/python

#################################################################################################################################################################################################################
#This script is an adaptation of Adaptor_Trimming_Cram2fastq_v0.5.py for general use.
#Lara Bossini-Castillo 13th October 2017. Wellcome Trust Sanger Institute
#lbc@sanger.ac.uk
#################################################################################################################################################################################################################

#PLEASE CHECK PATHS TO REFERENCE FILES AND SOFTWARE!!!!


USAGE='''python Adaptor_Trimming_Cram2fastq_vXX.py FILE.fastq --- This program is used for adapter triming and QC of ATAC-seq fastq files (paired-end sequencing), ;
usage:
    python %s [--prefix=#] [--ref] inputdata(fastq)
#input: fastq files input, it must obey the regular expression "\.((fastq)$"
}.log files
#defaults:
--prefix: Default: the basename of inputdata.
--ref: New alignment reference file
'''

# Import built-in modules
import os # using operating system dependent functionality
import os.path # common pathname manipulations
import getopt
import sys # access to some variables used or maintained by the interpreter and to functions that interact strongly with the interpreter
import time # time-related functions
import re # provides regular expression matching operations similar to those found in Perl
import subprocess #keeps command output

#Compression of fastq files
def gzip1(prefix,indir):
    newpath = r'input'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        print 'input directory created'
        cmd = 'gzip -f -c %s > %s' % (prefix+'.1.fastq',indir+'/input/'+prefix+'.1.fastq.gz')
        print cmd
        os.system(cmd)
    else:
        cmd = 'gzip -f -c %s > %s' % (prefix+'.1.fastq',indir+'/input/'+prefix+'.1.fastq.gz')
        print cmd
        os.system(cmd)
    return prefix+'.1.fastq'

def gzip2(prefix,indir):
    newpath = r'input'
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        print 'input directory created'
        cmd = 'gzip -f -c %s > %s' % (prefix+'.2.fastq',indir+'/input/'+prefix+'.2.fastq.gz')
        print cmd
        os.system(cmd)
    else:
        cmd = 'gzip -f -c %s > %s' % (prefix+'.2.fastq',indir+'/input/'+prefix+'.2.fastq.gz')
        print cmd
        os.system(cmd)
    return prefix+'.2.fastq'

#Removal of remaining Nextera adaptors

def adaptortrimming(gzip_1,gzip_2,ATAC_trim_adapters_directory, prefix, indir, primers_file, sampleIndexMap_file):
    newpath2 = r'output'
    if not os.path.exists(newpath2):
        os.makedirs(newpath2)
        print 'output folder created'
        cmd = 'echo %s | python %s/ATAC_trim_adapters_Adapted_v1.py --fastqDirIn %s --fastqDirOut %s --atacPrimers %s --sampleIndexMap %s' % (prefix, ATAC_trim_adapters_directory, indir+'/input/', indir+'/output/', primers_file, sampleIndexMap_file)
        print cmd
        os.system(cmd)
    else:
        cmd = 'echo %s | python %s/ATAC_trim_adapters_Adapted_v1.py --fastqDirIn %s --fastqDirOut %s --atacPrimers %s --sampleIndexMap %s' % (prefix, ATAC_trim_adapters_directory, indir+'/input/', indir+'/output/', primers_file, sampleIndexMap_file)
        print cmd
        os.system(cmd)
    return


#Realigment to GRCh38 build of the human genome
def realignment(adaptortrimming, bwa_directory, prefix, indir, ref_fa):
    cmd = '%s/bwa mem -M -t 12 %s %s %s > %s' % (bwa_directory, ref_fa, indir+'/output/'+prefix+'.1.trimmed.fastq.gz', indir+'/output/'+prefix+'.2.trimmed.fastq.gz', prefix+'.trimmed.ReAligned_GRCh38.sam')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sam'

#Output sam file to bam file
def bamCreation(realign,samtool_directory, prefix):
    cmd = '%s/samtools view -bS %s -o %s' % (samtool_directory, prefix+'.trimmed.ReAligned_GRCh38.sam', prefix+'.trimmed.ReAligned_GRCh38.bam')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.bam'

# Sorting of bam file
def bamsort(bamcreate,samtool_directory, prefix):
    cmd = '%s/samtools sort %s %s' % (samtool_directory, prefix+'.trimmed.ReAligned_GRCh38.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted.bam'

#Calls script adapted from previous work by Dafni Glinos to remove mitochondrial reads, duplicate reads and non-uniquely mapped reads
def ATAC_preQC(bamfile_realigned,atac_process_directory, prefix):
    cmd = 'python %s/atac_process_AdaptedGRCh38_v1.py %s --NOrm' % (atac_process_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted.bam')
    print cmd
    os.system(cmd)
    return

# Get read removal summary
def bamstats(repreQC,bamtool_directory,prefix):
    cmd = '%s/bamtools stats -in %s > %s' % (bamtool_directory,prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nomt_nodup.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nomt_nodup_bamstats.txt')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nomt_nodup_bamstats.txt'

#Bam file sort based on read name
def bamsort2(bamstat,samtool_directory, prefix):
    cmd = '%s/samtools sort -n %s %s' % (samtool_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nomt_nodup.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.bam'

#Bam file to bed file
def bam2bedpe(bamsort_2,bedtool_directory,prefix):
    cmd = '%s/bedtools bamtobed -bedpe -i %s > %s' % (bedtool_directory,prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.bedpe')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.bedpe'

#Call modified script based on Kaur Alasoo's https://github.com/kauralasoo/Blood_ATAC/blob/master/scripts/bedpe2bed.py to get insert fragment bed file (this step is needed for correct peak calling since BAMPE option in MACS2 does not work properly in ATAC-seq data. Specific shift and extension parameters need to be provided to MACS2 in this case.
#You can call peaks either at the fragment level or at the cut site level (individual transposase integration site). By including Kaur Alasoo's script https://github.com/kauralasoo/Blood_ATAC/blob/master/scripts/fragmentsToCutSites.py This option is not included by default.

def bedpe2bed(bam_2_bedpe, ATAC_trim_adapters_directory, prefix):
    cmd = 'less %s | python %s/bedpe2bed.py > %s' % (prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.bedpe', ATAC_trim_adapters_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.bed')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nomt_nodup.noPE.bed'


#Bed file sorting based on coordinates 
def sortbed(bedpe_2_bed, bedtool_directory, prefix):
    cmd = '%s/sortBed -i %s > %s' % (bedtool_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.bed', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.bed')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.bed'


#Peak calling using MACS2 and specific shift and extension parameters based on 'Alasoo et al. 2017. Genetic effects on chromatin accessibility foreshadow gene expression changes in macrophage immune response. BioRxiv'
def callPeak(sortBed,prefix):
    cmd = 'macs2 callpeak -t %s -f AUTO -n %s -g hs --nomodel --shift -25 --extsize 50 -B' % (prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.bed', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted')
    print cmd
    os.system(cmd)
    return

#Checking peak number
def peakCount(callPeaks, prefix):
    cmd = 'wc -l < %s' % (prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted_peaks.narrowPeak')
    output = subprocess.check_output(cmd, shell=True)
    if float(output) <= 10000:
        print prefix+' has less than 10,000 peaks'
        return
    else:
        print  prefix+' has '+output+' peaks.'
        return

#Cheching signal to noise ratio (also named FRiP). FRiP < 10% is considered low 
def signal2noiseRatio(peakCounts,bedtool_directory,prefix):
    cmd = '%s/bedtools intersect -a %s -b %s -c > %s' % (bedtool_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted_peaks.narrowPeak', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.bed', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.ReadsinPeaks')
    os.system(cmd)
    cmd2 = 'wc -l < %s' % (prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.bed')
    output = subprocess.check_output(cmd2, shell=True)
    cmd3 = 'awk \'{sum += $11} END {print sum}\' %s' % (prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nomt_nodup.noPE.sorted.ReadsinPeaks')
    output2 = subprocess.check_output(cmd3, shell=True)
    if float(output2) <= float(output)/10: #This division by 10 is due to a QC threshold of 10% FRiP that we established experimentaly
        print prefix+' has low signal to noise ratio. Number of reads in bed '+output+'. Number of reads in peaks '+output2
        exit()
    else:
        print  prefix+' has good signal to noise ratio. Number of reads in bed '+output+'. Number of reads in peaks '+output2
        return

def bamcount(bamfile):
    cmd = 'samtools view %s | wc -l ' % bamfile
    print cmd
    return ossystemresult(cmd)

def ossystemresult(command):
    fp=os.popen(command,'r') # Open a pipe to or from command. The return value is an open file object connected to the pipe
    return fp.read().rstrip()

if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 1: # The list of command line arguments passed to a Python script
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ['ref']) # arguments passed to the command line
    inputfile = args[0]

    regular = r'\.(fastq)$'
    if not re.search(regular, inputfile):
        print '''the input file must obey the regular expression "\.(fastq)$" '''
        sys.exit(1)
    postfix = re.search(regular, inputfile).group(0)
    
    #default:
    ref_fa = #Include your FILE.fa BWA reference genome
    prefix = re.sub(regular, "", os.path.basename(inputfile)) # Return the string obtained by replacing the leftmost non-overlapping occurrences of pattern in string by the replacement repl Return the base name of pathname path
    rmswitch = False
    
    bamtool_directory = #Include your bamtools path
    bedtool_directory = #Include your bedtools path
    ATAC_trim_adapters_directory = #Include path to adapter trimming script
    indir = os.getcwd()
    primers_file = 'ATAC_index_primers.txt' #File provided, add path if needed
    sampleIndexMap_file = 'Example_Sample_i5_i7_tags_3.txt' #Dummy file provided this is specific for you samples
    bwa_directory = #Include your BWA path
    samtool_directory = #Include your samtools path
    atac_process_directory = #Include path to ATAC-seq QC script
    genomefile_directory = 'hg38ChromSizes.genome' #Provided file obtained from https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
    
    for o, a in opts:
        if o == '--ref':
            ref_fa = a
        if o == '--prefix':
            prefix = a

    logfile = open(prefix+".log", "w")
    print >> logfile, 'inputfile\t%s' % (inputfile)

    gzip_1 = gzip1(prefix,indir)
    gzip_2 = gzip2(prefix,indir)
    adaptortrim = adaptortrimming(gzip_1,gzip_2,ATAC_trim_adapters_directory, prefix, indir, primers_file, sampleIndexMap_file)
    realign = realignment(adaptortrim,bwa_directory, prefix, indir, ref_fa)
    bamcreate = bamCreation(realign,samtool_directory, prefix)
    bamfile_realigned = bamsort(bamcreate,samtool_directory, prefix)
    repreQC = ATAC_preQC(bamfile_realigned,atac_process_directory, prefix)
    bamstat = bamstats(repreQC,bamtool_directory,prefix)
    bamsort_2 = bamsort2(bamstat,samtool_directory, prefix)
    bam_2_bedpe = bam2bedpe(bamsort_2,bedtool_directory,prefix)
    bedpe_2_bed = bedpe2bed(bam_2_bedpe, ATAC_trim_adapters_directory, prefix)
    sort_bed = sortbed(bedpe_2_bed, bedtool_directory, prefix)
    callingPeaks = callPeak(sort_bed,prefix)
    peakCounts = peakCount(callingPeaks, prefix)
    signal2noiseRatio(peakCounts,bedtool_directory,prefix)
    
    logfile.close()

    endtime = time.time()
    print 'it takes %.3f minutes or %.3f seconds' % ((endtime-starttime)/60, endtime -starttime)
