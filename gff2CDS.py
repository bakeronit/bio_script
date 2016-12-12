#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse,gzip
from collections import defaultdict
parser = argparse.ArgumentParser(description='Extract CDS from gff file and genome fasta file, can translate to protein',
                                 epilog="bakeronit@gmail.com")
parser.add_argument('gff', type=str, help='gff file')
parser.add_argument('genome',type=str, help='genome file')
parser.add_argument('-t','--translation',action='store_true',default = False,help='translate to protein or not,default:do not translate')
parser.add_argument('--longest',action='store_true',default = False,help='extract the longest cds or not,default do not output longest cds')
parser.add_argument('-p','--prefix',type=str,default = 'my',help='prefix of output file,default "my" ')
#TODO:Add some function
args = parser.parse_args()

def parseGff(file):
    genes = defaultdict(list)
    with open(file,'rt') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#') or line == '':
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            scaffoldid = cols[0]
            geneCDS = dict()
            if cols[2] == 'gene':
                geneid = cols[8].split(';')[0].split('=')[1]
                continue
            if cols[2].lower() == 'cds':
                mRNAid = cols[8].split(';')[1].split('=')[1]
                geneCDS['id'], geneCDS['strand'], geneCDS['start'],geneCDS['end'] = "%s|%s"%(geneid,mRNAid),\
                cols[6],int(cols[3]),int(cols[4])
                genes[scaffoldid].append(geneCDS)
    return genes

genes = parseGff(args.gff)

def reverseComplement(seq):
    complement = {'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T': 'A','t':'a','N':'N','n':'n'}
    reverse = ''
    for base in seq:
        reverse = complement[base] + reverse
    return reverse

CDS = defaultdict(str)
fh = gzip.open(args.genome,'rt') if args.genome.endswith('gz') else open(args.genome,'rt')
begin = 0
for line in fh:
    line = line.strip()
    if line.startswith('>'):
        if begin == 0:
            begin = 1
        elif not name in genes:
            pass
        else:
            for geneCDS in genes[name]:
                cdsId,strand,start,end = geneCDS['id'],geneCDS['strand'],geneCDS['start'],geneCDS['end']
                CDS[cdsId] = CDS[cdsId] + genome[start-1:end] if strand == '+' else CDS[cdsId] + reverseComplement(genome[start-1:end])
        name = line[1:]
        genome = ''
    else:
        genome += line
#last scaffold
for geneCDS in genes[name]:
    cdsId,strand,start,end = geneCDS['id'],geneCDS['strand'],geneCDS['start'],geneCDS['end']
    CDS[cdsId] = CDS[cdsId] + genome[start-1:end] if strand == '+' else CDS[cdsId] + reverseComplement(genome[start-1:end])
#delete big var
genes = genome = None
fh.close()
        
def translate(CDS,id):
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    pep = ''
    if not CDS.startswith('ATG'):
        print("CODON WARNING:%s not start with ATG"%id)
    if not len(CDS)%3 == 0:
        print("lENGTH WANRING:length of %s is not triple"%id)
    for n in range(0,len(CDS),3):
        if CDS[n:n+3] in codontable:
            pep += codontable[CDS[n:n+3]]
        else:
            pep += 'X'    
    return pep
             
out = open('%s.CDS.fa'%args.prefix,'wt')
out2 = open('%s.PEP.fa'%args.prefix,'wt')
for k in CDS:      
    print('>'+ k,file = out)
    if args.translation:
        print('>'+ k,file = out2)
        print(translate(CDS[k],k),file = out2)
    print(CDS[k], file = out)
out.close()
out2.close()

if args.longest:
    longestCDS = defaultdict(str)
    for k in CDS:
        kgeneid = k.split('|')[0]
        if len(CDS[k]) > len(longestCDS[kgeneid]):
            longestCDS[kgeneid] = CDS[k]

    out = open('%s.longestCDS.fa'%args.prefix,'wt')
    out2 = open('%s.longestPEP.fa'%args.prefix,'wt')
    for g in longestCDS:
        print('>' + g,file = out)
        if args.translation:
            print('>' + g,file = out2)
            print(translate(longestCDS[g],g),file = out2)
        print(longestCDS[g],file = out)
    out.close()
    out2.close()