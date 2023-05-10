#!/extdata6/Sukjun/opt/python3/bin/python3

#Script for Converting Old Sequence Data to Pseudo Fastq Files
#Data obtained from Barbiarz et. al. 2008 Genes Dev Paper
#January 9th 2015


bDEBUGMODE = 0

import os, sys, itertools, pickle, time, array, re, copy, subprocess, uuid

sBASE_DIR     = 'E:\Jinman\Programming\DecayElement\mESC'

nMIN_READ_LEN = 15
nMAX_READ_LEN = 50

DATASET1_2_FILENAME_SIZE = 3
DATASET1_2_READ_ID_SIZE  = 2
DATASET3_READLINE_SIZE   = 2


def main():
    for i in [1,2,3]:
        convert2fastq (i)
#def END: main


def convert2fastq (nDataset):
    print('Converting Dataset%s to Fastq %s' % (nDataset, time.ctime()))

    sDataDir  = '%s/Dataset%s' % (sBASE_DIR, nDataset if bDEBUGMODE==0 else 'Test')
    sOutDir   = '%s/fastq'     % sBASE_DIR
    os.makedirs(sOutDir, exist_ok=True)

    OutFile   = open('%s/SR0_Dataset%s.fastq' % (sOutDir, nDataset), 'w')

    if nDataset != 3: ## Dataset 1 or 2

        list_sFiles = os.listdir(sDataDir)
        nReadIDCnt  = 0

        for sFile in list_sFiles:
            if not sFile.endswith('.fa'): continue  #SKIP Files with invalid extension

            list_sFileName = sFile.split('.')

            if len(list_sFileName) < DATASET1_2_FILENAME_SIZE: continue #SKIP Files with invalid filename

            nReadLen       = int(list_sFileName[1])

            if not (nMIN_READ_LEN <= nReadLen <= nMAX_READ_LEN): continue #SKIP Fasta files of invalid read size


            sFileDir            = '%s/%s' % (sDataDir, sFile)
            OutFile, nReadIDCnt = parse_fa_file_dataset1_and_2 (sFileDir, OutFile, nReadIDCnt, nReadLen)

        #loop END: sFile

    else: # Dataset3: Differently formatted file than Dataset 1 and 2

        nReadIDCnt          = 0

        sFileDir            = '%s/GSM314552-10947.txt' % sDataDir

        OutFile, nReadIDCnt = parse_fa_file_dataset3 (sFileDir, OutFile, nReadIDCnt)

    #if END: nDataset

    print('Converting Dataset%s to Fastq %s...DONE' % (nDataset, time.ctime()))

    OutFile.close()
#def END: convert2fastq


def parse_fa_file_dataset1_and_2 (sFileDir, OutFile, nReadIDCnt, nReadLen):
    print('Parsing FA File %s %s' % (sFileDir, time.ctime()))

    InFile = open('%s' % (sFileDir), 'r')

    for sReadLine in InFile:
        #Fa File Format:
        #Line1:  ReadID   >4874_1
        #Line2:  ReadSeq  GGAATTGTGAGCGGA
        list_sReadID  = sReadLine.strip('\n').split('_')

        #V-S Check: List Size
        if len(list_sReadID) > DATASET1_2_READ_ID_SIZE:
            sys.exit('ERROR: list_sReadID Size= %d\nlist_sReadID= %s' % (len(list_sReadID), list_sReadID))

        sReadID     = list_sReadID[0]
        nReplicate  = int(list_sReadID[1])

        sReadSeq    = InFile.readline().strip('\n')

        #V-S Check: Read Length
        if len(sReadSeq) != nReadLen:
            sys.exit('ERROR: Invalid Read Size= %d File Designated Read Size=%d' % (len(sReadSeq), nReadLen))

        for i in range(nReplicate):
            nReadIDCnt += 1
            sNewReadID = '@%s_%s' % (sReadID[1:], nReadIDCnt)
            sReadQual  = 'I' * len(sReadSeq)

            OutFile.write('%s\n' % sNewReadID)
            OutFile.write('%s\n' % sReadSeq)
            OutFile.write('+\n')
            OutFile.write('%s\n' % sReadQual)
        #loop END: i
    #loop END: sReadLine
    InFile.close()

    print('Parsing FA File %s %s...DONE' % (sFileDir, time.ctime()))
    return OutFile, nReadIDCnt
#def END: parse_fa_file_dataset1_and_2


def parse_fa_file_dataset3 (sFileDir, OutFile, nReadIDCnt):
    print('Parsing Text File %s %s' % (sFileDir, time.ctime()))

    InFile = open(sFileDir, 'r')

    for sReadLine in InFile:
        #File Format:
        #Column:    1                       |   2
        # ID:       ReadSeq                 |   Replicate
        # Example:  ATAGAATATAACCTTTGCGTGT	|   1

        if sReadLine.startswith('#'): continue  #SKIP Comment Lines
        if sReadLine.startswith('ID'): continue #SKIP Column Labels

        list_sColumn        = sReadLine.strip('\n').split('\t')

        #V-S Check: Parsed Column Size
        if len(list_sColumn) > DATASET3_READLINE_SIZE:
            sys.exit('ERROR: list_sColumn Size= %d' % len(list_sColumn))

        sReadSeq, nReplicate = list_sColumn

        for i in range(int(nReplicate)):
            nReadIDCnt += 1

            sReadID   = '@DATA3_%s' % nReadIDCnt
            sReadQual = 'I' * len(sReadSeq)

            OutFile.write('%s\n' % sReadID)
            OutFile.write('%s\n' % sReadSeq)
            OutFile.write('+\n')
            OutFile.write('%s\n' % sReadQual)
        #loop END: i
    #loop END: sReadLine
    InFile.close()
    print('Parsing Text File %s %s...DONE' % (sFileDir, time.ctime()))

    return OutFile, nReadIDCnt

#def END: parse_fa_file_dataset3


main()
