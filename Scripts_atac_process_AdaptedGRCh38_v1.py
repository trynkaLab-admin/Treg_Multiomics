#!/usr/bin/python

USAGE = '''atac_process.py --- This program is used for aligned bam files from ATAC-seq;
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


# Remove mitochondrial reads
def bamMTrm(bamfile):
    cmd = 'samtools view -b %s chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chr1_KI270706v1_random chr1_KI270707v1_random chr1_KI270708v1_random chr1_KI270709v1_random chr1_KI270710v1_random chr1_KI270711v1_random chr1_KI270712v1_random chr1_KI270713v1_random chr1_KI270714v1_random chr2_KI270715v1_random chr2_KI270716v1_random chr3_GL000221v1_random chr4_GL000008v2_random chr5_GL000208v1_random chr9_KI270717v1_random chr9_KI270718v1_random chr9_KI270719v1_random chr9_KI270720v1_random chr11_KI270721v1_random chr14_GL000009v2_random chr14_GL000225v1_random chr14_KI270722v1_random chr14_GL000194v1_random chr14_KI270723v1_random chr14_KI270724v1_random chr14_KI270725v1_random chr14_KI270726v1_random chr15_KI270727v1_random chr16_KI270728v1_random chr17_GL000205v2_random chr17_KI270729v1_random chr17_KI270730v1_random chr22_KI270731v1_random chr22_KI270732v1_random chr22_KI270733v1_random chr22_KI270734v1_random chr22_KI270735v1_random chr22_KI270736v1_random chr22_KI270737v1_random chr22_KI270738v1_random chr22_KI270739v1_random chrY_KI270740v1_random chrUn_KI270302v1 chrUn_KI270304v1 chrUn_KI270303v1 chrUn_KI270305v1 chrUn_KI270322v1 chrUn_KI270320v1 chrUn_KI270310v1 chrUn_KI270316v1 chrUn_KI270315v1 chrUn_KI270312v1 chrUn_KI270311v1 chrUn_KI270317v1 chrUn_KI270412v1 chrUn_KI270411v1 chrUn_KI270414v1 chrUn_KI270419v1 chrUn_KI270418v1 chrUn_KI270420v1 chrUn_KI270424v1 chrUn_KI270417v1 chrUn_KI270422v1 chrUn_KI270423v1 chrUn_KI270425v1 chrUn_KI270429v1 chrUn_KI270442v1 chrUn_KI270466v1 chrUn_KI270465v1 chrUn_KI270467v1 chrUn_KI270435v1 chrUn_KI270438v1 chrUn_KI270468v1 chrUn_KI270510v1 chrUn_KI270509v1 chrUn_KI270518v1 chrUn_KI270508v1 chrUn_KI270516v1 chrUn_KI270512v1 chrUn_KI270519v1 chrUn_KI270522v1 chrUn_KI270511v1 chrUn_KI270515v1 chrUn_KI270507v1 chrUn_KI270517v1 chrUn_KI270529v1 chrUn_KI270528v1 chrUn_KI270530v1 chrUn_KI270539v1 chrUn_KI270538v1 chrUn_KI270544v1 chrUn_KI270548v1 chrUn_KI270583v1 chrUn_KI270587v1 chrUn_KI270580v1 chrUn_KI270581v1 chrUn_KI270579v1 chrUn_KI270589v1 chrUn_KI270590v1 chrUn_KI270584v1 chrUn_KI270582v1 chrUn_KI270588v1 chrUn_KI270593v1 chrUn_KI270591v1 chrUn_KI270330v1 chrUn_KI270329v1 chrUn_KI270334v1 chrUn_KI270333v1 chrUn_KI270335v1 chrUn_KI270338v1 chrUn_KI270340v1 chrUn_KI270336v1 chrUn_KI270337v1 chrUn_KI270363v1 chrUn_KI270364v1 chrUn_KI270362v1 chrUn_KI270366v1 chrUn_KI270378v1 chrUn_KI270379v1 chrUn_KI270389v1 chrUn_KI270390v1 chrUn_KI270387v1 chrUn_KI270395v1 chrUn_KI270396v1 chrUn_KI270388v1 chrUn_KI270394v1 chrUn_KI270386v1 chrUn_KI270391v1 chrUn_KI270383v1 chrUn_KI270393v1 chrUn_KI270384v1 chrUn_KI270392v1 chrUn_KI270381v1 chrUn_KI270385v1 chrUn_KI270382v1 chrUn_KI270376v1 chrUn_KI270374v1 chrUn_KI270372v1 chrUn_KI270373v1 chrUn_KI270375v1 chrUn_KI270371v1 chrUn_KI270448v1 chrUn_KI270521v1 chrUn_GL000195v1 chrUn_GL000219v1 chrUn_GL000220v1 chrUn_GL000224v1 chrUn_KI270741v1 chrUn_GL000226v1 chrUn_GL000213v1 chrUn_KI270743v1 chrUn_KI270744v1 chrUn_KI270745v1 chrUn_KI270746v1 chrUn_KI270747v1 chrUn_KI270748v1 chrUn_KI270749v1 chrUn_KI270750v1 chrUn_KI270751v1 chrUn_KI270752v1 chrUn_KI270753v1 chrUn_KI270754v1 chrUn_KI270755v1 chrUn_KI270756v1 chrUn_KI270757v1 chrUn_GL000214v1 chrUn_KI270742v1 chrUn_GL000216v2 chrUn_GL000218v1 chrEBV > %s' % (
    bamfile, bamfile.replace('.bam', '_nomt.bam'))
    os.system(cmd)
    return bamfile.replace('.bam', '_nomt.bam')


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
    uniqnoMTfile = bamMTrm(bamfile)
    uniqnoMTnodupbamfile = bam2nodupbam(uniqnoMTfile)
    bedfile = (uniqnoMTnodupbamfile, prefix)

    print >> logfile, 'input_reads\t%s' % bamcount(inputfile)
    print >> logfile, 'unique_reads\t%s' % bamcount(bamfile)
    print >> logfile, 'no_mt_reads\t%s' % bamcount(uniqnoMTfile)
    print >> logfile, 'nonredundant_reads\t%s' % bamcount(uniqnoMTnodupbamfile)
    logfile.close()

    if rmswitch:
        os.remove(bamfile)  # Remove (delete) the file path.
        os.remove(uniqnoMTfile)

    endtime = time.time()
    print 'it takes %.3f minutes or %.3f seconds' % ((endtime - starttime) / 60, endtime - starttime)
