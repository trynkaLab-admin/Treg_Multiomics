# Obtained from https://github.com/kauralasoo/Blood_ATAC/blob/master/scripts/bedpe2bed.py by Kaur Alasoo.
import argparse
import fileinput

parser = argparse.ArgumentParser(description="Convert BEDPE into a BED file of fragments.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--maxFragmentLength", help="Maximum fragment length between two pairs of reads.", default="1000")
args = parser.parse_args()

for line in fileinput.input("-"):
    line = line.rstrip()
    fields = line.split("\t")

    # Check for chimeras
    if fields[0] == fields[3]:
        if fields[0] != ".":
            fragment_length = int(fields[5]) - int(fields[1])
            if fragment_length < int(args.maxFragmentLength):
                out = "\t".join([fields[0], fields[1], fields[5], fields[6], str(fragment_length)])
                print(out)
