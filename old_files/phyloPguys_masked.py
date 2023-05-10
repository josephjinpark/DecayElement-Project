#! /extdata2/Sukjun/opt/python3/bin/python3
import sys,time,re,math,array,pickle
from collections import defaultdict
import itertools
CHRDIR = '/extdata3/Daekwan/data/Genome/Chromosomes/hg19/%s.fa'

def site_8m(sMiRNA) :
    return [site_7m8(sMiRNA)[0] + "A"]

def site_7m8(sMiRNA) :
    return [sMiRNA[-8:][:-1]]

def site_7a1(sMiRNA) :
    return [site_6m(sMiRNA)[0] + "A"]
def site_6m(sMiRNA) :
    return [sMiRNA[-7:][:-1]]

def fetch_cbc(in_file, start, end):
    in_f = open(in_file, 'rb')
    in_f.seek(start * 4)
    array_cbc = array.array('f', [])
    array_cbc.fromfile(in_f, end - start)
    in_f.close()
    #if strand == '-': array_cbc = array_cbc[::-1]
    return array_cbc

def fetch_cbc_cons(ref, scope, chrom, start, end):
    in_file_cbc = '/extdata2/Sukjun/project/biogenesis/cons/%s/%s/%s.cbc'%(ref, scope, chrom)
    cons = list(fetch_cbc(in_file_cbc, start, end))
    return cons

def wigfix2cbc(chromlen, in_file_wigfix, out_file_cbc):
    array_cbc = array.array('f', [0.0] * int(chromlen))
    for line in open(in_file_wigfix):
        if line.startswith('fixedStep'):
            id    = line.strip('\n').split()
            start = int(id[2].split('=')[1])  # 1-based
            step  = int(id[3].split('=')[1])
            assert step == 1
        else:
            array_cbc[start-1] = float(line.strip('\n'))
            start += 1
    array_cbc.tofile(open(out_file_cbc, 'wb'))


class RefSeq:
    def __init__(self):
        self.sGeneSym       = ""
        self.sNMID          = ""
        self.sChrID         = ""
        self.sStrand        = ""
        self.nTxnStart      = ""
        self.nTxnEnd        = ""
        self.nExonNum       = 0
        self.nORFStart      = 0
        self.nORFEnd        = 0
        self.nExonStartlist = []
        self.nExonEndlist   = []
        self.nExonSize      = 0
        self.sExonSeq       = ""
        self.n5UTRSize      = 0
        self.s5UTRSeq       = ""
        self.nORFSize       = 0
        self.sORFSeq        = ""
        self.n3UTRSize      = 0
        self.s3UTRSeq       = ""
        self.sRefSeqline    = ""
        self.fExonphyloP    = []
        self.f5UTRphyloP    = []
        self.fORFphyloP     = []
        self.f3UTRphyloP    = []

    def parse_ref_line(self,sRefline):
                #RefFlat string contents:InfoList[num]
                #0              1           2           3                   4                   
                #Gene symbol    RefSeqID    Chromosome  Strand direction    Txn start position
                #RTP3           NM_031440   chr3    +               46539484     
                #5                  6           7           8               9                       10
                #Txn end position   CDS start   CDS end     Number of exons Exon start positions    Exon end positions
                #46542439       46539552    46542389    2               46539484,46541845,      46539707,46542439,
        self.sRefSeqline    = sRefline.strip().replace(' ','\t')
        sInfoList           = self.sRefSeqline.split('\t')
        if len(sInfoList) != 11:
            print(self.sRefSeqline)
            sys.exit(1)
        
        self.sGeneSym       = sInfoList[0]
        self.sNMID          = sInfoList[1]
        self.sChrID         = sInfoList[2]  ##chr1~22,X,Y
        self.sStrand        = sInfoList[3]
        self.nTxnStart      = int(sInfoList[4])
        self.nTxnEnd        = int(sInfoList[5])
        self.nORFStart      = int(sInfoList[6])
        self.nORFEnd        = int(sInfoList[7])
        self.nExonNum       = int(sInfoList[8])
        self.nExonStartlist = [int(i) for i in sInfoList[9].strip(',').split(',')]
        self.nExonEndlist   = [int(i) for i in sInfoList[10].strip(',').split(',')]
        if (not(self.nExonNum == len(self.nExonStartlist) and self.nExonNum == len(self.nExonEndlist))) or (min(self.nExonStartlist+self.nExonEndlist) <= 0):
            print('Exonlist Not valid',self.sNMID,self.sGeneSym)
            sys.exit(1)
    ##parsing end

    def seqanalyze(self,sChrSeq):
        self.sExonSeq       = ''.join([sChrSeq[self.nExonStartlist[i]:self.nExonEndlist[i]] for i in range(self.nExonNum)])
        self.nExonSize      = len(self.sExonSeq)
        
        nORFStartExon       = 0
        while self.nORFStart > self.nExonEndlist[nORFStartExon]:
            self.n5UTRSize  += (self.nExonEndlist[nORFStartExon] - self.nExonStartlist[nORFStartExon]) 
            nORFStartExon   += 1
        self.n5UTRSize      += (self.nORFStart - self.nExonStartlist[nORFStartExon])

        nORFEndExon = self.nExonNum - 1
        while self.nORFEnd < self.nExonStartlist[nORFEndExon]:
            self.n3UTRSize  += (self.nExonEndlist[nORFEndExon] - self.nExonStartlist[nORFEndExon])
            nORFEndExon     -= 1
        self.n3UTRSize      += (self.nExonEndlist[nORFEndExon] - self.nORFEnd)

        self.nORFSize = self.nExonSize - self.n5UTRSize - self.n3UTRSize

        self.s5UTRSeq   = self.sExonSeq[:self.n5UTRSize]
        self.sORFSeq    = self.sExonSeq[self.n5UTRSize:(self.nExonSize-self.n3UTRSize)]
        self.s3UTRSeq   = self.sExonSeq[(self.nExonSize-self.n3UTRSize):]

        if self.sStrand == '-':
            self.sExonSeq   = complementaryReverse(self.sExonSeq)
            sTempUTR        = complementaryReverse(self.s3UTRSeq)
            self.s3UTRSeq   = complementaryReverse(self.s5UTRSeq)
            self.sORFSeq    = complementaryReverse(self.sORFSeq)
            self.s5UTRSeq   = sTempUTR
            self.n5UTRSize  = len(self.s5UTRSeq)
            self.n3UTRSize  = len(self.s3UTRSeq)
    ##seqanalyze end

    def phyloPassign(self,ref,scope):
        if self.sExonSeq == "":
            print('Check if seqanalyze module had been operated')
            sys.exit(1)
        fConslist_raw = fetch_cbc_cons(ref, scope, self.sChrID, self.nTxnStart, self.nTxnEnd)
        fConslist_frag= [fConslist_raw[(self.nExonStartlist[i]-self.nTxnStart):(self.nExonEndlist[i]-self.nTxnStart)] for i in range(self.nExonNum)]
        for frag in fConslist_frag:
            self.fExonphyloP = self.fExonphyloP + frag

        if self.sStrand =='-':
            self.fExonphyloP = self.fExonphyloP[::-1]
        
        self.f5UTRphyloP = self.fExonphyloP[:self.n5UTRSize]
        self.fORFphyloP  = self.fExonphyloP[self.n5UTRSize:(self.nExonSize-self.n3UTRSize)]
        self.f3UTRphyloP = self.fExonphyloP[(self.nExonSize-self.n3UTRSize):]
    ##phyloPassign end
    
    def canonical_mask(self, conservedlist,stype):
        temp_u3seq = self.s3UTRSeq
        typedic = {'8m':0,'7m8':1,'7A1':2,'6m':3,'6A1':4}

        for mirseq in conservedlist:
            mirseq_comrev = complementaryReverse(mirseq.upper().replace('U','T'))
            oc_match = site_8m(mirseq_comrev)[0]
            m8_match = site_7m8(mirseq_comrev)[0]
            A1_match = site_7a1(mirseq_comrev)[0]
            hx_match = site_6m(mirseq_comrev)[0]
            matchlist = [oc_match,m8_match,A1_match,hx_match]

            for match in matchlist[:typedic[stype]]:#[oc_match,m8_match,A1_match,hx_match]: #mask except 6m #[oc_match, m8_match, A1_match, hx_match]:
                for i in re.finditer(match,temp_u3seq):
                    maskseq = 'N'*len(match)
                    temp_u3seq = temp_u3seq[:i.start()]+maskseq+temp_u3seq[i.end():]
        self.s3UTRSeq = temp_u3seq
    ##canonical_mask end


def mean(values):
    if len(values) == 0:
        return None
    return sum(values, 0.0) / len(values)

def stderror(binlist):
    try:
        return standardDeviation(binlist,0.0)/math.sqrt(len(binlist))*1.96
    except TypeError:
        return 'Nan'

def standardDeviation(values, option):
    if len(values) < 2:
        return None

    sd = 0.0
    sum = 0.0
    meanValue = mean(values)

    for i in range(0, len(values)):
        diff = values[i] - meanValue
        sum += diff * diff

    sd = math.sqrt(sum / (len(values) - option))
    return sd


def complement(base):
    complbase = ''
    if base == 'A':
            complbase = 'T'
    elif base == 'T':
            complbase = 'A'
    elif base == 'G':
            complbase = 'C'
    elif base == 'C':
            complbase = 'G'
    return complbase

def complementaryReverse(sSeq):                                                                 ##Make complementary 5'-3' sequence
        sTempSeq = sSeq[::-1]
        retSeq = ''
        for base in sTempSeq:
                retSeq = retSeq + complement(base)
        return retSeq

def sequence(location):
        seqlist = open(location,"r")
        sTitle  = seqlist.readline()
        sSeq    = seqlist.read()
        if (sTitle[0]!=">"):
            sys.exit("Maybe this file is not a FASTA file." )
        else:
            sjoin = sSeq.replace("\n","").upper()
        seqlist.close()
        return sjoin

def genedic(refdir,stype):
    cRefDict = defaultdict(list)
    conservedlist = [i.strip().split('\t')[1] for i in open('./mir87_fullseq.txt','r')]
    with open(refdir,'r') as InF:
        reflines = [i.strip() for i in InF]
    for line in reflines:
        cGene = RefSeq()
        cGene.parse_ref_line(line)
        cRefDict[cGene.sChrID].append(cGene)
    for chrid in cRefDict:
        print(chrid,time.ctime())
        sChrSeq = sequence(CHRDIR %(chrid))
        for gene in cRefDict[chrid]:
            gene.seqanalyze(sChrSeq)
            gene.phyloPassign('hg19','vertebrate_100way')
            gene.canonical_mask(conservedlist,stype)
            #print(gene.sGeneSym)
            #print(mean(gene.f5UTRphyloP),mean(gene.fORFphyloP),mean(gene.f3UTRphyloP))
            #with open('./PLEKHS.txt' ,'w') as OutF:
            #    print(gene.s5UTRSeq, file = OutF)
            #    print(*gene.f5UTRphyloP, sep = '\t',file = OutF)
            #    print(gene.sORFSeq, file = OutF)
            #    print(*gene.fORFphyloP,sep = '\t',file = OutF)
            #    print(gene.s3UTRSeq, file = OutF)
            #    print(*gene.f3UTRphyloP, sep = '\t', file = OutF)
    with open('./ref/RefDic_100way.87masked.for%spickle' %stype,'wb') as OutF:
        pickle.dump(cRefDict,OutF)
for STYPE in ['8m','7m8','7A1','6m','6A1']:
    genedic('./ref/refFlat-12-Sep-2011_non_redundant.txt',STYPE)

#with open('./fafai.txt','r') as InF:
#    nums = [i.strip().split() for i in InF]
#    nums = [(i[0],int(i[1])) for i in nums]
#for num in nums:
#    print(num, time.ctime())
#    wigfix2cbc(num[1], '/extdata2/Sukjun/project/biogenesis/cons/hg19/vertebrate_100way/%s.phyloP100way.wigFix' %num[0],'/extdata2/Sukjun/project/biogenesis/cons/hg19/vertebrate_100way/%s.cbc' %num[0])
