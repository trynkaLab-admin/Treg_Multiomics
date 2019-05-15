#!/usr/bin/env python
import sys
import os

#PARAMETERS
executing_script = sys.argv[1]
print executing_script
filelist = sys.argv[2]
filelist_line = int(sys.argv[3]) - 1
f = open(filelist)
lines = f.readlines()
filelist_name = lines[filelist_line]
filelist_name = filelist_name.rstrip()
print filelist_name
command = "python " + executing_script + " " + filelist_name
print command
os.system(command)

