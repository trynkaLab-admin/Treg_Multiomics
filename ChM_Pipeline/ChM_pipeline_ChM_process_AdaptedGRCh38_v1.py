#!/usr/bin/python

USAGE = '''ChM_process_AdaptedGRCh38_v1.py --- This program is used for aligned bam files from ChM-seq;
usage:
    python %s [--prefix=#] [--NOrm] inputdata(bam)
#input: bam file, it must obey the regular expression "\.((bam)$"
#output: ${prefix}.bam, ${prefix}.bam.bai, ${prefix}.bed, ${prefix}.bw, ${prefix}.log files
#defaults:
#--ref: DNA reference chromosome sizes genome file.
--prefix: Default: the basename of inputdata.
--NOrm: if set the option, the program will not remove the bam file in the middle steps. Default: False.
'''
# Import built-in modules
import os  # using operating system dependent functionality
import os.path  # common pathname manipulations
import getopt
import sys  # access to some variables used or maintained by the interpreter and to functions that interact strongly
# with the interpreter
import time  # time-related functions
import re  # provides regular expression matching operations similar to those found in Perl


def ossystemresult(command):
    fp = os.popen(command,
                  'r')  # Open a pipe to or from command. The return value is an open file object connected to the pipe
    return fp.read().rstrip()


# Keep only uniquely mapped reads
def bam2uniqbam(bamfile, prefix):
    cmd = 'samtools view -bq 1 %s > %s' % (bamfile, prefix + '_uniq.bam')
    os.system(cmd)
    return prefix + '_uniq.bam'


# Remove duplicated reads
def bam2nodupbam(bamfile):
    cmd = 'samtools rmdup %s %s' % (bamfile, bamfile.replace('.bam', '_nodup.bam'))
    os.system(cmd)
    return bamfile.replace('.bam', '_nodup.bam')


# Create Bam index
def bamindex(bamfile):
    cmd = 'samtools index %s' % bamfile
    os.system(cmd)


# Counting reads in Bam file
def bamcount(bamfile):
    cmd = 'samtools view %s | wc -l ' % bamfile
    return ossystemresult(cmd)


def bamtobed(bamfile, prefix):
    cmd = 'bedtools bamtobed -i %s -split > %s' % (bamfile, prefix + '.bed')
    os.system(cmd)
    return prefix + '.bed'


if __name__ == '__main__':

    starttime = time.time()
    if len(sys.argv) < 1:  # The list of command line arguments passed to a Python script
        print USAGE % sys.argv[0]
        sys.exit(1)

    opts, args = getopt.getopt(sys.argv[1:], "", ['NOrm', 'ref', 'prefix='])  # arguments passed to the command line
    inputfile = args[0]

    regular = r'\.(bam)$'
    if not re.search(regular, inputfile):
        print '''the input file must obey the regular expression "\.(bam)$" '''
        sys.exit(1)
    postfix = re.search(regular, inputfile).group(0)

    # default: Return the string obtained by replacing the leftmost non-overlapping occurrences of pattern in string
    # by the replacement repl Return the base name of pathname path
    
    ref_fa = '/lustre/scratch117/cellgen/teamtrynka/lara/Temp_Resources/BWA_GRCh38_15_index/Homo_sapiens.GRCh38_15.fa'
    prefix = re.sub(regular, "", os.path.basename(inputfile))
    rmswitch = True

    for o, a in opts:
        if o == '--prefix':
            prefix = a
        if o == "--NOrm":
            rmswitch = False

    logfile = open(prefix + ".log", "w")
    print >> logfile, 'inputfile\t%s' % (inputfile)

    bamfile = bam2uniqbam(inputfile, prefix)
    bamindex(bamfile)
    uniqnoMTnodupbamfile = bam2nodupbam(bamfile)
    bedfile = (uniqnoMTnodupbamfile, prefix)

    print >> logfile, 'input_reads\t%s' % bamcount(inputfile)
    print >> logfile, 'unique_reads\t%s' % bamcount(bamfile)
    print >> logfile, 'nonredundant_reads\t%s' % bamcount(uniqnoMTnodupbamfile)
    logfile.close()

    if rmswitch:
        os.remove(bamfile)  # Remove (delete) the file path.

    endtime = time.time()
    print 'it takes %.3f minutes or %.3f seconds' % ((endtime - starttime) / 60, endtime - starttime)
