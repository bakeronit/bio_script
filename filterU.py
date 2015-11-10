#! /usr/bin/env python
from __future__ import print_function
import os
import argparse
import sys
parser= argparse.ArgumentParser(description='filter U in CDS.fasta')
parser.add_argument('-i','--input',type=str, help = 'input fasta - must given')
parser.add_argument('-l','--list',type=str,help='list file for input fasta')
parser.add_argument('-o','--outdir',type=str,default='./',help='outfile file directory, default is ./')
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args=parser.parse_args()

def readFasta(filename):
    U_file={}
    seq={}
    with open(filename,'rt') as fh:
        for line in fh:
            line=line.rstrip()
            if line.startswith('>'):
                ContainU = False
                seqname = line
                U_file[seqname] = 0
                seq[seqname]=list()
                continue
            else:
                for base in line:
                    if base == 'U' or base == 'u':
                        ContainU = True
                        U_file[seqname] += 1
                if ContainU == False:
                    seq[seqname].append(line)
                else:
                    if seqname in seq:
                        del seq[seqname]
                    pass
    return U_file,seq
if args.list and args.input:
    print("cannot have both -i and -l!")
    sys.exit()

if args.list:
    f=open(args.list,'rt')
    for fpath in f.readlines():
        fpath=fpath.rstrip()
        filename= os.path.basename(fpath)
        fh2=open(os.path.join(args.outdir,'TotalStat.out'),'a')
        U_file,seq=readFasta(fpath)
        total=0
        total_seq=0
        for i in U_file:
            if U_file[i]>0:
                total+=U_file[i]
                total_seq+=1
        sp_name=filename.split('.')[0]
        print('\nSpecie name:'+sp_name,file =fh2)
        print('>>Total seq contain U:%d'%total_seq, file = fh2)
        print('>>Total U:%d'%total , file= fh2)
        print('>>After filter:%d'%len(seq),file =fh2)
        if not total_seq==0:
            fh3=open(os.path.join(args.outdir,filename+'.filterd'),'wt')
            for i in seq:
                print(i,file= fh3)
                for s in seq[i]:
                    print(s,file = fh3)
if args.input:
    f=open(args.input,'rt')
    filename= os.path.basename(args.input)
    fh2=open(os.path.join(args.outdir,filename+'.stat.out'),'wt')
    U_file,seq=readFasta(args.input)
    total=0
    total_seq=0
    for i in U_file:
        if U_file[i]>0:
            total+=U_file[i]
            total_seq+=1
    print('>>Total seq contain U:%d'%total_seq, file = fh2)
    print('>>Total U:%d'%total , file= fh2)
    print('>>After filter:%d'%len(seq),file =fh2)
    filename= os.path.basename(args.input)
    if not total_seq==0:
        fh3=open(os.path.join(args.outdir,filename+'.filterd'),'wt')
        for i in seq:
            print(i,file= fh3)
            for s in seq[i]:
                print(s,file = fh3)



        
                    
