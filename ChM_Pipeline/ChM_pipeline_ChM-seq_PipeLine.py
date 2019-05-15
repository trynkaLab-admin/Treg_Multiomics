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

#Calls script adapted from previous work by Dafni Glinos to remove duplicate reads and non-uniquely mapped reads
def ATAC_preQC(bamfile_realigned,atac_process_directory, prefix):
    cmd = 'python %s/ChM_process_AdaptedGRCh38_v1.py %s --NOrm' % (atac_process_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted.bam')
    print cmd
    os.system(cmd)
    return

# Get read removal summary
#def bamstats(repreQC,bamtool_directory,prefix):
def bamstats(bamtool_directory,prefix):
    cmd = '%s/bamtools stats -in %s > %s' % (bamtool_directory,prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nodup.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nodup_bamstats.txt')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nodup_bamstats.txt'

#Bam file sort based on read name
def bamsort2(bamstat,samtool_directory, prefix):
    cmd = '%s/samtools sort -n %s %s' % (samtool_directory, prefix+'.trimmed.ReAligned_GRCh38.sorted_uniq_nodup.bam', prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup')
    print cmd
    os.system(cmd)
    return prefix+'.trimmed.ReAligned_GRCh38.sorted2_uniq_nodup.bam'

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
    ref_fa = '/lustre/scratch117/cellgen/teamtrynka/lara/Temp_Resources/BWA_GRCh38_15_index/Homo_sapiens.GRCh38_15.fa'
    prefix = re.sub(regular, "", os.path.basename(inputfile)) # Return the string obtained by replacing the leftmost non-overlapping occurrences of pattern in string by the replacement repl Return the base name of pathname path
    rmswitch = False
    
    bamtool_directory = '/software/svi/bin/bamtools-2.4.0/bin'
    bedtool_directory = '/software/hgi/pkglocal/bedtools-2.22.0/bin'
    ATAC_trim_adapters_directory = '/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM' 
    indir = os.getcwd()
    primers_file = '/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM/ChM_index_primers.txt'
    sampleIndexMap_file = 'ChM_sample_index_primers.txt'
    bwa_directory = '/software/hgi/pkglocal/bwa-0.7.9a/bin'
    samtool_directory = '/software/CGP/external-apps/samtools-0.1.9/bin'
    atac_process_directory = '/lustre/scratch117/cellgen/teamtrynka/lara/TregProject_ChM'
    genomefile_directory = 'hg38ChromSizes.genome'
    
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
    bamstat = bamstats(bamtool_directory,prefix)
    bamsort_2 = bamsort2(bamstat,samtool_directory, prefix)

    
    logfile.close()

    endtime = time.time()
    print 'it takes %.3f minutes or %.3f seconds' % ((endtime-starttime)/60, endtime -starttime)
