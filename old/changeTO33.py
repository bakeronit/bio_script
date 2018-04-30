#!/usr/bin/env python
# -*- coding: utf-8 -*-
#filename: change to phred33
from __future__ import print_function
import sys
import os
import argparse
import time

parser= argparse.ArgumentParser(description='change Phred+66 to Phred+33')
parser.add_argument('-i','--input',type=str, nargs='+',help = 'input fastq - must given')
parser.add_argument('-o','--outdir',type=str,default='./',help='outfile file directory, default is ./')
args=parser.parse_args()

def GETQUALITY(filename):
    with open(filename,'rt') as fh:
        row=0
        qmin=qmax=0
        for line in fh:
            line = line.rstrip()
            if row%4 == 3:
                for i in range(len(line)):
                    qmin = ord(line[i]) if qmin== 0 or ord(line[i]) < qmin else qmin
                    qmax = ord(line[i]) if ord(line[i]) > qmax else qmax
            row += 1
    return qmin,qmax

def Phred64toPhred33(filename,system,outdir):
    if system == 'Illumina (1.5+)':
        Q = 33
    elif system == 'Illumina (1.3+)':
        Q = 31
    elif system == 'Solexa':
        Q = 26
    elif system == 'Illumina (1.8+)' or 'Sanger':
        return 
    out=open(os.path.join(outdir,(os.path.basename(filename)+'.changed')),'wt')
    with open(filename,'rt') as fh:
        row=0
        for line in fh:
            line=line.rstrip()
            if row%4 == 3:
                change_line=''
                for i in range(len(line)):
                    change_line += chr(ord(line[i])-Q)
                    change_line.strip()
                print(change_line,file = out)
            else:
                 print(line,file = out)
            row += 1
    return 0
############### main function ###################
ISOTIMEFORMAT = '%Y-%m-%d %X'
print("start time:" + time.strftime(ISOTIMEFORMAT,time.localtime()))  
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
fh2=open(os.path.join(args.outdir,'Total_file.log'),'wt')
if os.path.isdir(args.input[0]):
    file_list=[os.path.join(args.input[0],name) for name in os.listdir(args.input[0])]
else:
    file_list=args.input
for f in file_list:
    qmin,qmax=GETQUALITY(f)
    if qmin >= 66 and qmin <= 104 and qmax >= 66 and qmax <= 104:
        system = 'Illumina (1.5+)'
        print(os.path.basename(f)+':Illumina (1.5+) (Phred+64, 66 to 104)',file = fh2)
    elif qmin >= 64 and qmin <= 104 and qmax >= 64 and qmax <= 104 :
        system = 'Illumina (1.3+)'
        print(os.path.basename(f)+':Illumina (1.3+) (Phred+64, 64 to 104)',file = fh2)
    elif qmin >= 59 and qmin <= 104 and qmax >= 59 and qmax <= 104 :
        system = 'Solexa'
        print(os.path.basename(f)+':Solexa (Phred+64, 59 to 104)',file = fh2)
    elif qmin >= 33 and qmin <= 73 and qmax >= 33 and qmax <= 73:
        system = 'Sanger'
        print(os.path.basename(f)+':Sanger (Phred+33, 33 to 73)',file = fh2)
    elif qmin >= 33 and qmin <= 74 and qmax >= 33 and qmax <= 74 :
        system = 'Illumina (1.8+)'
        print(os.path.basename(f)+':Illumina (1.8+) (Phred+33, 33 to 74)',file = fh2)
    Phred64toPhred33(f,system,args.outdir)
print("end time:" + time.strftime(ISOTIMEFORMAT,time.localtime()))           
