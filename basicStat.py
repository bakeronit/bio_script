#! usr/bin/env python
# -*- coding: utf-8 -*-
# filename:basicStat.py
import time
import argparse
import collections
import matplotlib
import os
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

parser= argparse.ArgumentParser(description='some basic statistics of fastq file')
parser.add_argument('input',type=str, help = 'input fastq - must given')
parser.add_argument('-o','--outdir',type=str,default='./',help='outfile file directory, default is ./')
args=parser.parse_args()
ISOTIMEFORMAT = '%Y-%m-%d %X'

def ASCIItoquality33(ch):
    return ord(ch)-33
#def ASCIItoquality64(ch):
 #   return ord(cd)-64

def basicStat(filename):
    baseCount = collections.Counter()   # count how many ATCGN or other abnormal in sequece
    with open (filename,'r') as f:
        row=0
        baseNum=0
        read_length=0
        for line in f:
            line = line.rstrip()
            if  row%4 == 1:
                if read_length == 0:
                    n=len(line)     #min length of reads
                    m=len(line)     #max length of reads
                    #read_length=len(line)
                    read_length=[n,m]
                    #gc = [0] * (read_length+20)   # get stastistic of gc content cross read at every base postion.  std ~ 20bp
                    #gc_len=[0] * (read_length + 20)
                    gc=[[0]*2 for l in range(read_length[1]+20)]             #gc[count_gc,count_base]
                    quality=[0]*60           
                baseNum+=len(line)
                baseCount.update(line)  # update counter for every read
                read_length[0]=len(line) if len(line) < read_length[0] else read_length[0]
                read_length[1]=len(line) if len(line) > read_length[1] else read_length[1]
                for i in range(len(line)):
                    gc[i][1]+=1
                    if line[i] == 'G' or line[i] == 'C':
                        gc[i][0]+=1
            if row%4 == 3:                 # quanlity stastistic
                for q in line:
                    quality[ASCIItoquality33(q)]+=1      # default is phred33
            row+=1   # after perform +1 , row equal the number of lines have been read.
    GC_content = (baseCount['G'] + baseCount['C'])*100/baseNum
    fh = open(os.path.join(args.outdir,'Stat.out'),'w')              #output file
    print('>>Total sequece:%d'%(row/4),file = fh)
    print('>>Read length:%d-%d'%(read_length[0],read_length[1]), file=fh)
    print('>>Total Base:%d'%baseNum , file = fh)
    print('>>Readable:%.2fGb\t%.2fMb'%((baseNum/10e9),(baseNum/10e6)), file = fh)
    print('>>Base Composition:', file = fh)
    for i in baseCount:
        print(i+":%d"%baseCount[i], file = fh)
    print('>>GC content:%.2f'%GC_content + '%', file = fh)
    #fh2 = open(args.outfile,'w+')
    gc_per_base=[0]*read_length[1]
    print('\n>>GC content across per base:', file = fh)
    for i in range(read_length[1]):
        gc_per_base[i] = gc[i][0]/gc[i][1]   # total line divide 4 equal number of reads.
        print('%d\t%.2f'%((i+1),(gc_per_base[i]*100))+'%', file = fh)
    return gc_per_base,read_length, filename,quality
  
########## main function #################
print("start time:" + time.strftime(ISOTIMEFORMAT,time.localtime()))
gc,read_length ,filename,quality= basicStat(args.input)
#print(quality)
import numpy as np
import matplotlib.pyplot as plt
plt.figure()
plt.plot(gc)
plt.xlim(0,read_length[1])
plt.ylim(0.2,0.8)
plt.xlabel('Per base')
plt.ylabel('GC%')
plt.title('GC content across all bases in '+ filename)
plt.savefig(os.path.join(args.outdir,'GC.png'))
plt.figure()
x=np.arange(len(quality))
plt.bar(x,quality,alpha= .5, color='g')
plt.savefig(os.path.join(args.outdir,'hist.png'))
n,bins,patches=plt.hist(quality,)
print("end time:" + time.strftime(ISOTIMEFORMAT,time.localtime()))
