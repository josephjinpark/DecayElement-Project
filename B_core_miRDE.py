#!/extdata6/Sukjun/opt/python3/bin/python3

bDEBUGMODE = 0

import os, sys, itertools, pickle, time, array, re, copy, subprocess, uuid

from collections import Counter
from scipy import stats


## region Globals

sBASE_DIR         = '/extdata5/Jinman/05_DecayElement'
sCONSERVATION_DIR = '/extdata5/Jinman/conservation/'
sRNAfold_DIR      = '/home/lab/bin/ViennaRNA-2.0.5/Progs/RNAsubopt'
sFILEID_DIR_sRNA  = '/extdata5/Jinman/00_CG_srSeq/05_File_IDs'
sFILEID_DIR_mRNA  = '/extdata5/Jinman/03_CG_mrSeq/05_fileids'


fPSEUDO_COUNT     = float('1.0e-300')

list_sCancerTypes = ['breast', 'kidney', 'prostate', 'liver', 'thyroid',
                     'lung_adeno', 'lung_squa', 'headneck', 'stomach',
                     'uterine', 'kidney_pap', 'kidney_chrom', 'bladder']

list_sSpecies     = ['hsa', 'mmu', 'mml', 'bta', 'rno', 'gal']

list_sCHRID_HSA   = [   '1','2','3','4','5','6','7','8','9','10',
                        '11','12','13','14','15','16','17','18','19','20',
                        '21', '22','X','Y']

dict_sGENOME_DIR  = {   'hsa':'/extdata5/Jinman/cancermir/genome/hg19/hg19.fa',
                        'mmu':'/extdata5/Jinman/reference_genome/mm9/mm9.fa',
                        'mml':'/extdata5/Jinman/reference_genome/MMUL1.0/MMUL1.0.fa',
                        'rno':'/extdata5/Jinman/reference_genome/rn4/rn4.fa',
                        'bta':'/extdata5/Jinman/reference_genome/bosTau6/bosTau6.fa',
                        'gal':'/extdata5/Jinman/reference_genome/galGal4/galGal4.fa',
                        'TCGA':'/extdata5/Jinman/cancermir/genome/hg19/hg19.fa'}

dict_sMIRBASE_DIR = {   'hsa': '/extdata5/Jinman/reference_genome/mirbase/hsa.liftover.gff3',
                        'mmu': '/extdata5/Jinman/reference_genome/mirbase/mmu.gff3',
                        'mml': '/extdata5/Jinman/reference_genome/mirbase/mml.liftover.gff3',
                        'bta': '/extdata5/Jinman/reference_genome/mirbase/bta.liftover.gff3',
                        'rno': '/extdata5/Jinman/reference_genome/mirbase/rno.liftover.gff3',
                        'gal': '/extdata5/Jinman/reference_genome/mirbase/gga.liftover.gff3'}

dict_nSEED_LENGTH  = {'8mer': 8, '7m8': 7, '7A1': 7, '6mer': 6}

nCONSERVED_MIRNAS             = 197

# Validity Check for Parsing Input Text
nGFF3_COLUMN_SIZE             = 9
nCONSERVED_MIRNA_COLUMN_SIZE  = 2

## Conservation Score Filters
nFIRST_CONS_FILTER            = 7
nSECOND_CONS_FILTER           = 20
nMAX_INTERDISTANCE            = 20

# 3' Additional Matching Size
nADD_MATCH_MER_SIZE           = 5
nMIN_ADDITIONAL_SITE_SIZE     = 6
nMAX_ADDITIONAL_SITE_SIZE     = 25

# Extension Length
nEXTEND_LEN                   = 50


sTIME_STAMP                   = '%s' % (time.ctime().replace(' ', '-').replace(':', '_') )

## endregion


## region Naming Convention:
#       n-      = int
#       s-      = string
#       f-      = float
#       c-      = class object
#       dict_-  = dictionaries
#       list_-  = lists
## endregion


## region Util Functions

## region SVG Plot Util Function
def svg2png(sInDir):
    list_sScript = []

    list_sScript.append('/extdata5/Jinman/bin/jre1.7.0/bin/java ')
    list_sScript.append('-jar /extdata5/Jinman/bin/batik-1.7/batik-rasterizer.jar ')
    list_sScript.append('%s > /dev/null' % sInDir)

    list_sScript = ''.join(list_sScript)

    os.system(list_sScript)
#def END: svg2png


def rnasubopt_e(seq, deltaEnergy=4.0):
    script = 'echo -ne "%s" | /home/lab/bin/ViennaRNA-2.0.5/Progs/RNAsubopt -s -e %d'%(seq, deltaEnergy)
    stdout = subprocess.Popen(script, stdout=subprocess.PIPE, shell=True).stdout
    list_dotbr = []
    for i, line in enumerate(stdout):
        line = str(line, 'UTF-8').strip('\n')
        if i >= 1:
            col = line.split()
            dotbr  = col[0]
            energy = float(col[1])
            list_dotbr.append([dotbr, energy])
    return list_dotbr
def adjust_coor(coor):
    coor_x = [ c[0] for c in coor ]
    coor_y = [ c[1] for c in coor ]
    min_x = min(coor_x)
    min_y = min(coor_y)
    adjcoor = [ (c[0]-min_x, c[1]-min_y) for c in coor ]
    return adjcoor
def load_ps(sInFile):

    in_f_ps = open(sInFile)
    fCoor_ori  = []
    nPairs_ori = []
    while True:
        sLine = in_f_ps.readline()
        if sLine.startswith('/coor ['):
            while True:
                sLine = in_f_ps.readline()
                if sLine.startswith('] def'): break
                fCoorX, fCoorY = map(lambda x: float(x), sLine[1:-2].split())
                fCoor_ori.append( (fCoorX, fCoorY) )
        elif sLine.startswith('/pairs ['):
            while True:
                sLine = in_f_ps.readline()
                if sLine.startswith('] def'): break
                nIdx1, nIdx2 = map(lambda x: int(x)-1, sLine[1:-2].split())
                nPairs_ori.append( (nIdx1, nIdx2) )
        elif sLine.startswith('%%EOF'): break
        else: pass
    in_f_ps.close()
    fCoor  = adjust_coor(fCoor_ori)
    nPairs = nPairs_ori
    return fCoor, nPairs
def rnaplot(sSeq, sDotbr):
    sUuid          = uuid.uuid4()
    sOut           = '>%s\n%s\n%s'  % (sUuid, sSeq, sDotbr)
    sTempDir       = '%s/temp'      % sBASE_DIR
    os.makedirs(sTempDir, exist_ok=True)

    sOutFile       = '%s.out'       % sUuid
    sOutFileDir    = '%s/%s'        % (sTempDir, sOutFile)
    sOutFile_ps    = '%s/%s_ss.ps'  % (sTempDir, sUuid)

    OutFile        = open(sOutFileDir, 'w')
    OutFile.write(sOut)
    OutFile.close()

    os.system('cd %s; /home/lab/bin/ViennaRNA-2.0.5/Progs/RNAplot < %s > /dev/null' % (sTempDir, sOutFile))

    fCoor, nPairs = load_ps(sOutFile_ps)
    os.system('rm -rf %s %s' % (sOutFileDir, sOutFile_ps))
    return fCoor
#def END: rnaplot

def strcheight(dotbr):
    h   = 0
    ret = []
    for d in dotbr:
        if   d == '(': h += 1; ret.append(h)
        elif d == ')': ret.append(h); h -= 1
        else: ret.append(h)
    return ret
def pairidx(dotbr, form=1):
    height   = strcheight(dotbr)
    list_idh = [ (i, d, h) for i, d, h in zip(range(len(dotbr)),dotbr,height) ]
    ret1     = [None for d in dotbr]
    ret2     = []
    for mh in range(max(height)+1):
        list_i = [ i for i, d, h in list_idh if h==mh and d=='(' ]
        list_j = [ i for i, d, h in list_idh if h==mh and d==')' ]
        for i, j in zip(list_i, list_j):
            ret1[i] = j
            ret1[j] = i
            ret2.append((i, j))
    return [ret1,ret2][form-1]

## endregion


## region Fetch Conservation Score
def fetch_cbc(sInDir, nStartPos, nEndPos):

    InFile    = open(sInDir, 'rb')

    InFile.seek(nStartPos * 4)

    array_cbc = array.array('f', [])

    array_cbc.fromfile(InFile, nEndPos - nStartPos)

    InFile.close()

    return array_cbc
#def fetch_cbc

def fetch_cbc_cons(sGenome, sScope, sChrID, nStartPos, nEndPos):

    InFileDir      = '%s/%s/%s/%s.cbc' % (sCONSERVATION_DIR, sGenome, sScope, sChrID)

    list_sConScore = list(fetch_cbc(InFileDir, nStartPos, nEndPos))

    return list_sConScore
#def END: fetch_cbc_cons


## endregion


## region Access ALN File

def get_overlap_length(nReadPos, nMiRNAPos):

    nStartPos_read, nEndPos_read = nReadPos
    nStartPos_mi, nEndPos_mi     = nMiRNAPos

    # V-S Check
    assert(nStartPos_read < nEndPos_read) and(nStartPos_mi < nEndPos_mi)

    # miRNA   :--------S--------------E-------
    # read    :------------S--------------E---
    # overlap :------------************-------

    nStartPos_overlap = max([nStartPos_read, nStartPos_mi])
    nEndPos_overlap   = min([nEndPos_read, nEndPos_mi])

    nOverlapLen = nEndPos_overlap - nStartPos_overlap + 1

    return nOverlapLen
#def END: get_overlap_length


def load_ali_file(sInputFile_ali):
    dict_sReadData = {}

    InFile         = open(sInputFile_ali, 'r')

    for sReadLine in InFile:
        list_sColumn = sReadLine.strip('\n').split('\t')

        sKey           = list_sColumn[0]
        nStartPos_read = int(list_sColumn[1])
        nEndPos_read   = int(list_sColumn[2])
        nPointer       = int(list_sColumn[4])

        if sKey not in dict_sReadData:
            dict_sReadData[sKey] = []

        dict_sReadData[sKey].append([nStartPos_read, nEndPos_read, nPointer])
    #loop END: sReadLine

    return dict_sReadData
#def END: load_ali_file


def access_aln_file(sInputFile_aln, sChrID, sStrand, nStartPos_miRNA, nEndPos_miRNA, nNorm, sRunType):

    dict_sReadData = load_ali_file(sInputFile_aln[:-4] + '.ali')
    sTargetKey    = sChrID + sStrand

    # Obtains the list of pointer values for the given chromosome and strand combo if there is overlap
    try: list_nPointer = [nPointer for nStartPos_read, nEndPos_read, nPointer in dict_sReadData[sTargetKey]
                         if get_overlap_length([nStartPos_read, nEndPos_read], [nStartPos_miRNA, nEndPos_miRNA]) > 0]
    except KeyError:
        list_nPointer = []
    #try END: list_nPointer

    nCoe          = pickle.load(open(sInputFile_aln[:-4]+'.coe','rb')) if nNorm else 1.0
    InFile        = open(sInputFile_aln, 'rb')

    if sRunType == 'expr':
        list_sReadData = generate_readlist_for_5p_expr (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA)
        return list_sReadData

    elif sRunType == 'muta':
        list_sReadData =  generate_readlist_for_variant_call (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA)
        return list_sReadData

    elif sRunType == 'reads':
        list_sReadData =  generate_readlist_for_reads_analysis (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA)
        list_sReadData = sorted(list_sReadData, key = lambda index : index[1])

        if sStrand == '-':
            list_sReadData = [ [reverse_complement(sSeq), nStartPos, nEndPos, nCnt, sType] for sSeq, nStartPos, nEndPos, nCnt, sType in list_sReadData]

        list_sReadData     = [['%s,%s,%s,%s' % (sSeq, nStartPos, nEndPos, sType), nCnt*nCoe] for sSeq, nStartPos, nEndPos, nCnt, sType in list_sReadData]
        dict_nReadCNT_RPM  = {}
        for sKey, fRPM in list_sReadData:
            if sKey not in dict_nReadCNT_RPM:
                dict_nReadCNT_RPM[sKey] = [0,0]

            dict_nReadCNT_RPM[sKey][0] += 1
            dict_nReadCNT_RPM[sKey][1] += fRPM
        #loop END: sKey, nCnt
        return dict_nReadCNT_RPM
    #if END: sRunType
    InFile.close()
#def END: access_aln_file


def generate_readlist_for_5p_expr (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA):

    list_sOutput = []
    for nPointer in list_nPointer:

        InFile.seek(nPointer)
        list_sReadData = pickle.load(InFile)

        for nStartPos_read, nEndPos_read, nReadCount, list_sUniq, list_sMult in list_sReadData:
            nUniqCnt       = 0
            nMultiCnt      = 0
            if get_overlap_length([nStartPos_read, nEndPos_read], [nStartPos_miRNA, nEndPos_miRNA]) > (nEndPos_read - nStartPos_read + 1) / 2:

                if sStrand == '+':
                    n5p_pos = nStartPos_read - nStartPos_miRNA
                    n3p_pos = nEndPos_read   - nStartPos_miRNA
                elif sStrand == '-':
                    n5p_pos = nEndPos_miRNA  - nEndPos_read
                    n3p_pos = nEndPos_miRNA  - nStartPos_read
                #if END: sStrand

                nUniqCnt  = sum([(nCnt * nCoe) for nCnt, sSeq in list_sUniq])
                nMultiCnt = sum([(nCnt * nCoe) for nCnt, sSeq in list_sMult])

                list_sOutput.append([n5p_pos, n3p_pos, nReadCount, nUniqCnt, nMultiCnt])
            #if END: get_overlap_length
        #loop END: nStartPos_read, nEndPos_read, nReadCount, nUniqCnt, nMultiCnt, sUniqSeqList, sMultSeqList
    #loop END: nPointer
    return list_sOutput
#def END: generate_readlist_for_expr


def generate_readlist_for_variant_call (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA):

    list_sOutput = []
    for nPointer in list_nPointer:

        InFile.seek(nPointer)
        list_sReadData = pickle.load(InFile)

        for nStartPos_read, nEndPos_read, nReadCount, list_sUniq, list_sMult in list_sReadData:
            if get_overlap_length([nStartPos_read, nEndPos_read], [nStartPos_miRNA, nEndPos_miRNA]) > (nEndPos_read - nStartPos_read + 1) / 2:

                list_sUniq = [[(nCnt*nCoe), sSeq] for nCnt, sSeq in list_sUniq]
                list_sMult = [[(nCnt*nCoe), sSeq] for nCnt, sSeq in list_sMult]

                list_sOutput.append([nStartPos_read, nEndPos_read, nReadCount, list_sUniq, list_sMult])
            #if END: get_overlap_length
        #loop END: nStartPos_read, nEndPos_read, nReadCount, nUniqCnt, nMultiCnt, sUniqSeqList, sMultSeqList
    #loop END: nPointer
    return list_sOutput
#def END: generate_readlist_for_variant_call


def generate_readlist_for_reads_analysis (InFile, list_nPointer, nCoe, sStrand, nStartPos_miRNA, nEndPos_miRNA):

    list_sOutput = []
    for nPointer in list_nPointer:

        InFile.seek(nPointer)
        list_sReadData = pickle.load(InFile)

        for nStartPos_read, nEndPos_read, nReadCount, list_sUniq, list_sMult in list_sReadData:
            if get_overlap_length([nStartPos_read, nEndPos_read], [nStartPos_miRNA, nEndPos_miRNA]) > (nEndPos_read - nStartPos_read + 1) / 2:

                list_sUniq = [[sSeq, int(nStartPos_read), int(nEndPos_read), nCnt, 'Uniq'] for nCnt, sSeq in list_sUniq]
                list_sMult = [[sSeq, int(nStartPos_read), int(nEndPos_read), nCnt, 'Mult'] for nCnt, sSeq in list_sMult]

                list_sOutput += list_sUniq + list_sMult
            #if END: get_overlap_length
        #loop END: nStartPos_read, nEndPos_read, nReadCount, nUniqCnt, nMultiCnt, sUniqSeqList, sMultSeqList
    #loop END: nPointer

    return list_sOutput
#def END: generate_readlist_for_reads_analysis

## endregion


def get_file_id_list(sFileName, bSeqType):
    if bSeqType == 'sRNA': sInDir = sFILEID_DIR_sRNA
    else: sInDir = sFILEID_DIR_mRNA
    return [id[:-1] for id in open('%s/fileid_%s.txt' %(sInDir, sFileName), 'r')]
#def END: get_file_id_list


def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    list_sSeq       = list(sSeq.upper()) # Turns the sequence in to a gigantic list
    list_sSeq       = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1] # Empty string, join the list of bases, [start : end : backwards]
#def END: reverse_complement


def include(sChrID_mir, sChrID_read, list_nArmPos, list_nMirPos):
    # coordinate is closed at both ends
    nStartPos_read, nEndPos_read = list_nArmPos
    nStartPos_mi, nEndPos_mi     = list_nMirPos

    #print('read', nStartPos_read, nEndPos_read)
    #print('mi', nStartPos_mi, nEndPos_mi)

    # V-S Check
    assert(nStartPos_read < nEndPos_read) and(nStartPos_mi < nEndPos_mi)
    if sChrID_mir.lower() != sChrID_read.lower():
        bFlag = 1
    else:

        if (nStartPos_mi <= nStartPos_read) and (nEndPos_mi >= nEndPos_read):
            bFlag = 0

        else:
            bFlag = 1
    #if END: sChrID

    return bFlag
#def END: include


def assign_cell_line_info (sFileName):
    dict_cell_line = {}
    ## Mayr Cell Line - sRNA-Seq
    dict_cell_line['GSM416733'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'HEK293T'}
    dict_cell_line['GSM416753'] = {'genome' : 'hsa', 'tissue_source': 'cervix', 'cell_line':'Hela'}
    dict_cell_line['GSM416754'] = {'genome' : 'hsa', 'tissue_source': 'bone', 'cell_line':'U2OS'}
    dict_cell_line['GSM416755'] = {'genome' : 'hsa', 'tissue_source': 'bone', 'cell_line':'143B'}
    dict_cell_line['GSM416756'] = {'genome' : 'hsa', 'tissue_source': 'lung', 'cell_line':'A549'}
    dict_cell_line['GSM416757'] = {'genome' : 'hsa', 'tissue_source': 'lung', 'cell_line':'H520'}
    dict_cell_line['GSM416758'] = {'genome' : 'hsa', 'tissue_source': 'colon', 'cell_line':'SW480'}
    dict_cell_line['GSM416759'] = {'genome' : 'hsa', 'tissue_source': 'colon', 'cell_line':'DLD2'}
    dict_cell_line['GSM416760'] = {'genome' : 'hsa', 'tissue_source': 'breast', 'cell_line':'MCF7'}
    dict_cell_line['GSM416761'] = {'genome' : 'hsa', 'tissue_source': 'breast', 'cell_line':'MB-MDA231'}

    ## ENCODE Cell Line  - sRNA-Seq
    dict_cell_line['SRR446301'] = {'genome' : 'hsa', 'tissue_source': 'cervix', 'cell_line':'Hela-S3'}
    dict_cell_line['SRR446302'] = {'genome' : 'hsa', 'tissue_source': 'cervix', 'cell_line':'Hela-S3'}
    dict_cell_line['SRR446303'] = {'genome' : 'hsa', 'tissue_source': 'brain', 'cell_line':'SK-N-SH_RA'}
    dict_cell_line['SRR446304'] = {'genome' : 'hsa', 'tissue_source': 'brain', 'cell_line':'SK-N-SH_RA'}
    dict_cell_line['SRR446306'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'H1-hESC'}
    dict_cell_line['SRR446307'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'H1-hESC'}
    dict_cell_line['SRR446308'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'H1-hESC'}
    dict_cell_line['SRR446309'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'AGO4450'}
    dict_cell_line['SRR446310'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'AGO4450'}
    dict_cell_line['SRR446315'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'BJ'}
    dict_cell_line['SRR446316'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'BJ'}
    dict_cell_line['SRR446317'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'NHEK'}
    dict_cell_line['SRR446318'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'NHEK'}
    dict_cell_line['SRR446319'] = {'genome' : 'hsa', 'tissue_source': 'cervix', 'cell_line':'Hela-S3'}
    dict_cell_line['SRR446320'] = {'genome' : 'hsa', 'tissue_source': 'cervix', 'cell_line':'Hela-S3'}
    dict_cell_line['SRR446323'] = {'genome' : 'hsa', 'tissue_source': 'breast', 'cell_line':'MCF7'}
    dict_cell_line['SRR446324'] = {'genome' : 'hsa', 'tissue_source': 'breast', 'cell_line':'MCF7'}
    dict_cell_line['SRR446329'] = {'genome' : 'hsa', 'tissue_source': 'lung', 'cell_line':'A549'}
    dict_cell_line['SRR446330'] = {'genome' : 'hsa', 'tissue_source': 'lung', 'cell_line':'A549'}
    dict_cell_line['SRR446337'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'NHEK'}
    dict_cell_line['SRR446338'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'NHEK'}
    dict_cell_line['SRR446388'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562'}
    dict_cell_line['SRR446389'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562'}
    dict_cell_line['SRR446392'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'GM12878'}
    dict_cell_line['SRR446393'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'GM12878'}
    dict_cell_line['SRR446394'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'GM12878'}
    dict_cell_line['SRR446395'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'GM12878'}
    dict_cell_line['SRR446398'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562'}
    dict_cell_line['SRR446399'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562'}

    ## BaekLab - Human Cell Line 12  - sRNA-Seq
    dict_cell_line['SR000A_001'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'Molt-4'}
    dict_cell_line['SR000A_002'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'THP-1'}
    dict_cell_line['SR000A_003'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'RPMI8226'}
    dict_cell_line['SR000A_004'] = {'genome' : 'hsa', 'tissue_source': 'muscle', 'cell_line':'A-673'}
    dict_cell_line['SR000A_005'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'ChangLiver'}
    dict_cell_line['SR000A_006'] = {'genome' : 'hsa', 'tissue_source': 'liver', 'cell_line':'Hep3B'}
    dict_cell_line['SR000A_007'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'HEK293T'}
    dict_cell_line['SR000A_008'] = {'genome' : 'hsa', 'tissue_source': 'kidney', 'cell_line':'A-704'}
    dict_cell_line['SR000A_009'] = {'genome' : 'hsa', 'tissue_source': 'normal', 'cell_line':'HEL299'}
    dict_cell_line['SR000A_010'] = {'genome' : 'hsa', 'tissue_source': 'lung', 'cell_line':'A549'}
    dict_cell_line['SR000A_011'] = {'genome' : 'hsa', 'tissue_source': 'breast', 'cell_line':'BT-20'}
    dict_cell_line['SR000A_012'] = {'genome' : 'hsa', 'tissue_source': 'testes', 'cell_line':'Tera-1'}

    ## BaekLab - Mouse Cell Line 12  - sRNA-Seq
    dict_cell_line['SR000A_013'] = {'genome' : 'mmu', 'tissue_source': 'lymphocyte', 'cell_line':'EL4'}
    dict_cell_line['SR000A_014'] = {'genome' : 'mmu', 'tissue_source': 'macrophage ', 'cell_line':'RAW254.7'}
    dict_cell_line['SR000A_015'] = {'genome' : 'mmu', 'tissue_source': 'lymphocyte', 'cell_line':'NS-1'}
    dict_cell_line['SR000A_016'] = {'genome' : 'mmu', 'tissue_source': 'fibroblast', 'cell_line':'NOR-10'}
    dict_cell_line['SR000A_017'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'NCTC'}
    dict_cell_line['SR000A_018'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'Hepa-1c1c7'}
    dict_cell_line['SR000A_019'] = {'genome' : 'mmu', 'tissue_source': 'kidney ', 'cell_line':'TCMK-1'}
    dict_cell_line['SR000A_020'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'Renca'}
    dict_cell_line['SR000A_021'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'MLg'}
    dict_cell_line['SR000A_022'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'LA-4'}
    dict_cell_line['SR000A_023'] = {'genome' : 'mmu', 'tissue_source': 'breast', 'cell_line':'MTV/TM-011'}
    dict_cell_line['SR000A_024'] = {'genome' : 'mmu', 'tissue_source': 'testes', 'cell_line':'F9'}

    ## BaekLab - Mouse Tissue 15  - sRNA-Seq
    dict_cell_line['SR0008_041'] = {'genome' : 'mmu', 'tissue_source': 'wholebrain', 'cell_line':'wholebrain'}
    dict_cell_line['SR0008_042'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon'}
    dict_cell_line['SR0008_043'] = {'genome' : 'mmu', 'tissue_source': 'heart', 'cell_line':'heart'}
    dict_cell_line['SR0008_044'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney'}
    dict_cell_line['SR0008_045'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver'}
    dict_cell_line['SR0008_046'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'lung'}
    dict_cell_line['SR0008_047'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary'}
    dict_cell_line['SR0008_048'] = {'genome' : 'mmu', 'tissue_source': 'pancreas', 'cell_line':'pancreas'}
    dict_cell_line['SR0008_049'] = {'genome' : 'mmu', 'tissue_source': 'skeletalmuscle', 'cell_line':'skeletalmuscle'}
    dict_cell_line['SR0008_050'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen'}
    dict_cell_line['SR0008_051'] = {'genome' : 'mmu', 'tissue_source': 'E10', 'cell_line':'E10'}
    dict_cell_line['SR0008_052'] = {'genome' : 'mmu', 'tissue_source': 'E12', 'cell_line':'E12'}
    dict_cell_line['SR0008_053'] = {'genome' : 'mmu', 'tissue_source': 'E13', 'cell_line':'E13'}
    dict_cell_line['SR0008_054'] = {'genome' : 'mmu', 'tissue_source': 'E15', 'cell_line':'E15'}
    dict_cell_line['SR0008_055'] = {'genome' : 'mmu', 'tissue_source': 'E17', 'cell_line':'E17'}

    ## GEO Data - Human Cell Line - mRNA-Seq
    dict_cell_line['SRR307903'] = {'genome' : 'hsa', 'tissue_source': 'skin', 'cell_line':'BJ_GEO'}
    dict_cell_line['SRR307904'] = {'genome' : 'hsa', 'tissue_source': 'skin', 'cell_line':'BJ_GEO'}
    dict_cell_line['SRR307897'] = {'genome' : 'hsa', 'tissue_source': 'lymph', 'cell_line':'GM12878_GEO'}
    dict_cell_line['SRR307898'] = {'genome' : 'hsa', 'tissue_source': 'lymph', 'cell_line':'GM12878_GEO'}
    dict_cell_line['SRR307911'] = {'genome' : 'hsa', 'tissue_source': 'ES', 'cell_line':'H1-hESC_GEO'}
    dict_cell_line['SRR307912'] = {'genome' : 'hsa', 'tissue_source': 'ES', 'cell_line':'H1-hESC_GEO'}
    dict_cell_line['SRR315336'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562_GEO'}
    dict_cell_line['SRR315337'] = {'genome' : 'hsa', 'tissue_source': 'blood', 'cell_line':'K562_GEO'}
    dict_cell_line['SRR315327'] = {'genome' : 'hsa', 'tissue_source': 'skin', 'cell_line':'NHEK_GEO'}
    dict_cell_line['SRR315328'] = {'genome' : 'hsa', 'tissue_source': 'skin', 'cell_line':'NHEK_GEO'}
    dict_cell_line['SRR315329'] = {'genome' : 'hsa', 'tissue_source': 'skin', 'cell_line':'NHEK_GEO'}
    dict_cell_line['SRR315315'] = {'genome' : 'hsa', 'tissue_source': 'brain', 'cell_line':'SK-N-SH_RA_GEO'}
    dict_cell_line['SRR315316'] = {'genome' : 'hsa', 'tissue_source': 'brain', 'cell_line':'SK-N-SH_RA_GEO'}

    ## GEO Data - Mouse Tissue 42 - mRNA-Seq
    dict_cell_line['SRR453166'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453167'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453168'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453169'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453170'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453171'] = {'genome' : 'mmu', 'tissue_source': 'colon', 'cell_line':'colon_GEO'}
    dict_cell_line['SRR453172'] = {'genome' : 'mmu', 'tissue_source': 'heart', 'cell_line':'heart_GEO'}
    dict_cell_line['SRR453173'] = {'genome' : 'mmu', 'tissue_source': 'heart', 'cell_line':'heart_GEO'}
    dict_cell_line['SRR453174'] = {'genome' : 'mmu', 'tissue_source': 'heart', 'cell_line':'heart_GEO'}
    dict_cell_line['SRR453175'] = {'genome' : 'mmu', 'tissue_source': 'heart', 'cell_line':'heart_GEO'}
    dict_cell_line['SRR453144'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453145'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453146'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453147'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453148'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453149'] = {'genome' : 'mmu', 'tissue_source': 'kidney', 'cell_line':'kidney_GEO'}
    dict_cell_line['SRR453150'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453151'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453152'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453153'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453154'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453155'] = {'genome' : 'mmu', 'tissue_source': 'liver', 'cell_line':'liver_GEO'}
    dict_cell_line['SRR453156'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'lung_GEO'}
    dict_cell_line['SRR453157'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'lung_GEO'}
    dict_cell_line['SRR453158'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'lung_GEO'}
    dict_cell_line['SRR453159'] = {'genome' : 'mmu', 'tissue_source': 'lung', 'cell_line':'lung_GEO'}
    dict_cell_line['SRR453077'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453078'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453079'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453080'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453081'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453082'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453083'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453084'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453085'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453086'] = {'genome' : 'mmu', 'tissue_source': 'ovary', 'cell_line':'ovary_GEO'}
    dict_cell_line['SRR453160'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}
    dict_cell_line['SRR453161'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}
    dict_cell_line['SRR453162'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}
    dict_cell_line['SRR453163'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}
    dict_cell_line['SRR453164'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}
    dict_cell_line['SRR453165'] = {'genome' : 'mmu', 'tissue_source': 'spleen', 'cell_line':'spleen_GEO'}

    ## GEO Data - Mouse Embryo - Sandberg Study - 40 Samples
    dict_cell_line['SRR805449'] = {'genome' : 'mmu', 'tissue_source': 'Zygote1', 'cell_line':'Zygote_1'}
    dict_cell_line['SRR805450'] = {'genome' : 'mmu', 'tissue_source': 'Zygote2', 'cell_line':'Zygote_3'}
    dict_cell_line['SRR805451'] = {'genome' : 'mmu', 'tissue_source': 'Zygote3', 'cell_line':'Zygote_3'}
    dict_cell_line['SRR805452'] = {'genome' : 'mmu', 'tissue_source': 'Zygote4', 'cell_line':'Zygote_4'}

    dict_cell_line['SRR805286'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_0-1', 'cell_line':'Early2cell_1'}
    dict_cell_line['SRR805287'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_0-1', 'cell_line':'Early2cell_2'}
    dict_cell_line['SRR805288'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_1-1', 'cell_line':'Early2cell_3'}
    dict_cell_line['SRR805289'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_1-2', 'cell_line':'Early2cell_4'}
    dict_cell_line['SRR805290'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_2-1', 'cell_line':'Early2cell_5'}
    dict_cell_line['SRR805291'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_2-2', 'cell_line':'Early2cell_6'}
    dict_cell_line['SRR805292'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_3-1', 'cell_line':'Early2cell_7'}
    dict_cell_line['SRR805293'] = {'genome' : 'mmu', 'tissue_source': 'Early2cell_3-2', 'cell_line':'Early2cell_8'}

    dict_cell_line['SRR805377'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_0-1', 'cell_line':'Mid2cell_1'}
    dict_cell_line['SRR805378'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_0-2', 'cell_line':'Mid2cell_2'}
    dict_cell_line['SRR805379'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_3-1', 'cell_line':'Mid2cell_3'}
    dict_cell_line['SRR805380'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_3-2', 'cell_line':'Mid2cell_4'}
    dict_cell_line['SRR805381'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_4-1', 'cell_line':'Mid2cell_5'}
    dict_cell_line['SRR805382'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_4-2', 'cell_line':'Mid2cell_6'}
    dict_cell_line['SRR805383'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_5-1', 'cell_line':'Mid2cell_7'}
    dict_cell_line['SRR805384'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_5-2', 'cell_line':'Mid2cell_8'}
    dict_cell_line['SRR805385'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_6-1', 'cell_line':'Mid2cell_9'}
    dict_cell_line['SRR805386'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_6-2', 'cell_line':'Mid2cell_10'}
    dict_cell_line['SRR805387'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_7-1', 'cell_line':'Mid2cell_11'}
    dict_cell_line['SRR805388'] = {'genome' : 'mmu', 'tissue_source': 'Mid2cell_8-2', 'cell_line':'Mid2cell_12'}

    dict_cell_line['SRR805337'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_5-1', 'cell_line':'Late2cell_1'}
    dict_cell_line['SRR805338'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_5-2', 'cell_line':'Late2cell_2'}
    dict_cell_line['SRR805339'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_6-1', 'cell_line':'Late2cell_3'}
    dict_cell_line['SRR805340'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_6-2', 'cell_line':'Late2cell_4'}
    dict_cell_line['SRR805341'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_7-1', 'cell_line':'Late2cell_5'}
    dict_cell_line['SRR805342'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_7-2', 'cell_line':'Late2cell_6'}
    dict_cell_line['SRR805343'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_8-1', 'cell_line':'Late2cell_7'}
    dict_cell_line['SRR805344'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_8-2', 'cell_line':'Late2cell_8'}
    dict_cell_line['SRR805345'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_9-1', 'cell_line':'Late2cell_9'}
    dict_cell_line['SRR805346'] = {'genome' : 'mmu', 'tissue_source': 'Late2cell_9-2', 'cell_line':'Late2cell_10'}

    dict_cell_line['SRR805294'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-1', 'cell_line':'Earlyblast_1'}
    dict_cell_line['SRR805295'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-10', 'cell_line':'Earlyblast_2'}
    dict_cell_line['SRR805296'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-12', 'cell_line':'Earlyblast_3'}
    dict_cell_line['SRR805297'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-15', 'cell_line':'Earlyblast_4'}
    dict_cell_line['SRR805298'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-16', 'cell_line':'Earlyblast_5'}
    dict_cell_line['SRR805299'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-17', 'cell_line':'Earlyblast_6'}
    dict_cell_line['SRR805300'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-18', 'cell_line':'Earlyblast_7'}
    dict_cell_line['SRR805301'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-2', 'cell_line':'Earlyblast_8'}
    dict_cell_line['SRR805302'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-22', 'cell_line':'Earlyblast_9'}
    dict_cell_line['SRR805303'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-3', 'cell_line':'Earlyblast_10'}
    dict_cell_line['SRR805304'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-4', 'cell_line':'Earlyblast_11'}
    dict_cell_line['SRR805305'] = {'genome' : 'mmu', 'tissue_source': 'Earlyblast_2-5', 'cell_line':'Earlyblast_12'}

    dict_cell_line['SRR805389'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-1', 'cell_line':'Midblast_1'}
    dict_cell_line['SRR805390'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-10', 'cell_line':'Midblast_2'}
    dict_cell_line['SRR805391'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-11', 'cell_line':'Midblast_3'}
    dict_cell_line['SRR805392'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-12', 'cell_line':'Midblast_4'}
    dict_cell_line['SRR805393'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-13', 'cell_line':'Midblast_5'}
    dict_cell_line['SRR805394'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-14', 'cell_line':'Midblast_6'}
    dict_cell_line['SRR805395'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-15', 'cell_line':'Midblast_7'}
    dict_cell_line['SRR805396'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-16', 'cell_line':'Midblast_8'}
    dict_cell_line['SRR805397'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-17', 'cell_line':'Midblast_9'}
    dict_cell_line['SRR805398'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-18', 'cell_line':'Midblast_10'}
    dict_cell_line['SRR805399'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-19', 'cell_line':'Midblast_11'}
    dict_cell_line['SRR805400'] = {'genome' : 'mmu', 'tissue_source': 'Midblast_1-2', 'cell_line':'Midblast_12'}

    dict_cell_line['SRR805347'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-10', 'cell_line':'Lateblast_1'}
    dict_cell_line['SRR805348'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-11', 'cell_line':'Lateblast_2'}
    dict_cell_line['SRR805349'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-12', 'cell_line':'Lateblast_3'}
    dict_cell_line['SRR805350'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-13', 'cell_line':'Lateblast_4'}
    dict_cell_line['SRR805351'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-16', 'cell_line':'Lateblast_5'}
    dict_cell_line['SRR805352'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-19', 'cell_line':'Lateblast_6'}
    dict_cell_line['SRR805353'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-2', 'cell_line':'Lateblast_7'}
    dict_cell_line['SRR805354'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-20', 'cell_line':'Lateblast_8'}
    dict_cell_line['SRR805355'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-21', 'cell_line':'Lateblast_9'}
    dict_cell_line['SRR805356'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-23', 'cell_line':'Lateblast_10'}
    dict_cell_line['SRR805357'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-24', 'cell_line':'Lateblast_11'}
    dict_cell_line['SRR805358'] = {'genome' : 'mmu', 'tissue_source': 'Lateblast_1-26', 'cell_line':'Lateblast_12'}

    dict_cell_line['SRR805223'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-1', 'cell_line':'4cell_1'}
    dict_cell_line['SRR805224'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-2', 'cell_line':'4cell_2'}
    dict_cell_line['SRR805225'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-4', 'cell_line':'4cell_3'}
    dict_cell_line['SRR805226'] = {'genome' : 'mmu', 'tissue_source': '4cell_2-1', 'cell_line':'4cell_4'}
    dict_cell_line['SRR805227'] = {'genome' : 'mmu', 'tissue_source': '4cell_2-1', 'cell_line':'4cell_5'}
    dict_cell_line['SRR805228'] = {'genome' : 'mmu', 'tissue_source': '4cell_2-3', 'cell_line':'4cell_6'}
    dict_cell_line['SRR805229'] = {'genome' : 'mmu', 'tissue_source': '4cell_2-4', 'cell_line':'4cell_7'}
    dict_cell_line['SRR805230'] = {'genome' : 'mmu', 'tissue_source': '4cell_3-1', 'cell_line':'4cell_8'}
    dict_cell_line['SRR805231'] = {'genome' : 'mmu', 'tissue_source': '4cell_3-3', 'cell_line':'4cell_9'}
    dict_cell_line['SRR805232'] = {'genome' : 'mmu', 'tissue_source': '4cell_3-4', 'cell_line':'4cell_10'}
    dict_cell_line['SRR805233'] = {'genome' : 'mmu', 'tissue_source': '4cell_4-1', 'cell_line':'4cell_11'}
    dict_cell_line['SRR805234'] = {'genome' : 'mmu', 'tissue_source': '4cell_4-2', 'cell_line':'4cell_12'}

    dict_cell_line['SRR805237'] = {'genome' : 'mmu', 'tissue_source': '8cell_1-1', 'cell_line':'8cell_1'}
    dict_cell_line['SRR805238'] = {'genome' : 'mmu', 'tissue_source': '8cell_1-2', 'cell_line':'8cell_2'}
    dict_cell_line['SRR805244'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-1', 'cell_line':'8cell_3'}
    dict_cell_line['SRR805245'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-2', 'cell_line':'8cell_4'}
    dict_cell_line['SRR805246'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-3', 'cell_line':'8cell_5'}
    dict_cell_line['SRR805247'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-4', 'cell_line':'8cell_6'}
    dict_cell_line['SRR805248'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-6', 'cell_line':'8cell_7'}
    dict_cell_line['SRR805249'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-7', 'cell_line':'8cell_8'}
    dict_cell_line['SRR805250'] = {'genome' : 'mmu', 'tissue_source': '8cell_2-8', 'cell_line':'8cell_9'}
    dict_cell_line['SRR805251'] = {'genome' : 'mmu', 'tissue_source': '8cell_5-1', 'cell_line':'8cell_10'}
    dict_cell_line['SRR805252'] = {'genome' : 'mmu', 'tissue_source': '8cell_5-2', 'cell_line':'8cell_11'}
    dict_cell_line['SRR805253'] = {'genome' : 'mmu', 'tissue_source': '8cell_5-3', 'cell_line':'8cell_12'}

    dict_cell_line['SRR805173'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-10', 'cell_line':'16cell_1'}
    dict_cell_line['SRR805174'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-11', 'cell_line':'16cell_2'}
    dict_cell_line['SRR805175'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-12', 'cell_line':'16cell_3'}
    dict_cell_line['SRR805176'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-13', 'cell_line':'16cell_4'}
    dict_cell_line['SRR805177'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-14', 'cell_line':'16cell_5'}
    dict_cell_line['SRR805178'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-15', 'cell_line':'16cell_6'}
    dict_cell_line['SRR805179'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-2', 'cell_line':'16cell_7'}
    dict_cell_line['SRR805180'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-3', 'cell_line':'16cell_8'}
    dict_cell_line['SRR805181'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-4', 'cell_line':'16cell_9'}
    dict_cell_line['SRR805182'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-5', 'cell_line':'16cell_10'}
    dict_cell_line['SRR805183'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-6', 'cell_line':'16cell_11'}
    dict_cell_line['SRR805184'] = {'genome' : 'mmu', 'tissue_source': '16cell_1-7', 'cell_line':'16cell_12'}

    ## GEO Data - Mouse Embryo - Zhong Study - 19 Samples
    dict_cell_line['SRR1267943'] = {'genome' : 'mmu', 'tissue_source': 'Zygote1', 'cell_line':'Zygote_1'}
    dict_cell_line['SRR1267944'] = {'genome' : 'mmu', 'tissue_source': 'Zygote2', 'cell_line':'Zygote_2'}
    dict_cell_line['SRR1267945'] = {'genome' : 'mmu', 'tissue_source': 'Zygote3', 'cell_line':'Zygote_3'}
    dict_cell_line['SRR1267946'] = {'genome' : 'mmu', 'tissue_source': 'Zygote4', 'cell_line':'Zygote_4'}

    dict_cell_line['SRR1267952'] = {'genome' : 'mmu', 'tissue_source': '2cell_1-1', 'cell_line':'2cell_1'}
    dict_cell_line['SRR1267953'] = {'genome' : 'mmu', 'tissue_source': '2cell_1-2', 'cell_line':'2cell_2'}
    dict_cell_line['SRR1267954'] = {'genome' : 'mmu', 'tissue_source': '2cell_2-1', 'cell_line':'2cell_3'}
    dict_cell_line['SRR1267955'] = {'genome' : 'mmu', 'tissue_source': '2cell_2-2', 'cell_line':'2cell_4'}

    dict_cell_line['SRR1267972'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-1', 'cell_line':'4cell_1'}
    dict_cell_line['SRR1267973'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-2', 'cell_line':'4cell_2'}
    dict_cell_line['SRR1267974'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-3', 'cell_line':'4cell_3'}
    dict_cell_line['SRR1267975'] = {'genome' : 'mmu', 'tissue_source': '4cell_1-4', 'cell_line':'4cell_4'}

    dict_cell_line['SRR1267992'] = {'genome' : 'mmu', 'tissue_source': 'InnerCell_1', 'cell_line':'InnerCell_1'}
    dict_cell_line['SRR1267993'] = {'genome' : 'mmu', 'tissue_source': 'InnerCell_2', 'cell_line':'InnerCell_2'}
    dict_cell_line['SRR1267994'] = {'genome' : 'mmu', 'tissue_source': 'InnerCell_3', 'cell_line':'InnerCell_3'}
    dict_cell_line['SRR1267995'] = {'genome' : 'mmu', 'tissue_source': 'InnerCell_4', 'cell_line':'InnerCell_4'}

    dict_cell_line['SRR1267996'] = {'genome' : 'mmu', 'tissue_source': 'Trophectoderm_1', 'cell_line':'Trophectoderm_1'}
    dict_cell_line['SRR1267997'] = {'genome' : 'mmu', 'tissue_source': 'Trophectoderm_2', 'cell_line':'Trophectoderm_2'}
    dict_cell_line['SRR1267998'] = {'genome' : 'mmu', 'tissue_source': 'Trophectoderm_3', 'cell_line':'Trophectoderm_3'}


    ## Barbiarz Paper - Mouse ES 3  - sRNA-Seq
    dict_cell_line['SR0_Dataset1'] = {'genome' : 'mmu', 'tissue_source': 'ES_454', 'cell_line':'ES_454'}
    dict_cell_line['SR0_Dataset2'] = {'genome' : 'mmu', 'tissue_source': 's_5_TagCount', 'cell_line':'s_5_TagCount'}
    dict_cell_line['SR0_Dataset3'] = {'genome' : 'mmu', 'tissue_source': 'ES_Illumina', 'cell_line':'ES_Illumina'}


    ## BaekLab - 19 Extended  - sRNA-Seq
    dict_cell_line['SR0004_001'] = {'genome' : 'hsa', 'tissue_source': 'cervix',   'cell_line':'Hela'}
    dict_cell_line['SR0004_002'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'CCRF-CEM'}
    dict_cell_line['SR0004_003'] = {'genome' : 'hsa', 'tissue_source': 'colon',    'cell_line':'HCT116'}
    dict_cell_line['SR0004_004'] = {'genome' : 'hsa', 'tissue_source': 'liver',    'cell_line':'HepG2'}
    dict_cell_line['SR0004_005'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'HL60'}
    dict_cell_line['SR0004_006'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'HuT78'}
    dict_cell_line['SR0004_007'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'IM9'}
    dict_cell_line['SR0004_008'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'Jiyoye'}
    dict_cell_line['SR0004_009'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'Jurkat'}
    dict_cell_line['SR0004_010'] = {'genome' : 'hsa', 'tissue_source': 'bone',     'cell_line':'KG1'}
    dict_cell_line['SR0004_011'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'K562'}
    dict_cell_line['SR0004_012'] = {'genome' : 'hsa', 'tissue_source': 'prostate', 'cell_line':'LNCaP'}
    dict_cell_line['SR0004_013'] = {'genome' : 'hsa', 'tissue_source': 'breast',   'cell_line':'MCF7'}
    dict_cell_line['SR0004_014'] = {'genome' : 'hsa', 'tissue_source': 'breast',   'cell_line':'MDA-MB-231'}
    dict_cell_line['SR0004_015'] = {'genome' : 'hsa', 'tissue_source': 'prostate', 'cell_line':'PC3'}
    dict_cell_line['SR0004_016'] = {'genome' : 'hsa', 'tissue_source': 'blood',    'cell_line':'Ramos'}
    dict_cell_line['SR0004_017'] = {'genome' : 'hsa', 'tissue_source': 'brain',    'cell_line':'SH-SY5Y'}
    dict_cell_line['SR0004_018'] = {'genome' : 'hsa', 'tissue_source': 'colon',    'cell_line':'SW620'}
    dict_cell_line['SR0004_019'] = {'genome' : 'hsa', 'tissue_source': 'brain',    'cell_line':'U87MG'}

    ## BaekLab - mRNA 14 Human - mRNA-Seq
    dict_cell_line['H42'] = {'genome' : 'hsa', 'tissue_source': 'cervix',   'cell_line':'Hela'}
    dict_cell_line['H43'] = {'genome' : 'hsa', 'tissue_source': 'normal',   'cell_line':'HEK293T'}
    dict_cell_line['H44'] = {'genome' : 'hsa', 'tissue_source': 'muscle',   'cell_line':'A-673'}
    dict_cell_line['H45'] = {'genome' : 'hsa', 'tissue_source': 'liver',    'cell_line':'Hep3B'}
    dict_cell_line['H46'] = {'genome' : 'hsa', 'tissue_source': 'prostate', 'cell_line':'PC3'}
    dict_cell_line['H47'] = {'genome' : 'hsa', 'tissue_source': 'colon',    'cell_line':'HCT116'}
    dict_cell_line['H48'] = {'genome' : 'hsa', 'tissue_source': 'breast',   'cell_line':'MDA-MB-231'}
    dict_cell_line['H49'] = {'genome' : 'hsa', 'tissue_source': 'breast',   'cell_line':'MCF7'}
    dict_cell_line['H50'] = {'genome' : 'hsa', 'tissue_source': 'lung',     'cell_line':'A549'}
    dict_cell_line['H51'] = {'genome' : 'hsa', 'tissue_source': 'prostate', 'cell_line':'LNCaP'}
    dict_cell_line['H52'] = {'genome' : 'hsa', 'tissue_source': 'colon',    'cell_line':'SW620'}
    dict_cell_line['H54'] = {'genome' : 'hsa', 'tissue_source': 'breast',   'cell_line':'U87-MG'}
    dict_cell_line['H57'] = {'genome' : 'hsa', 'tissue_source': 'liver',    'cell_line':'HepG2'}
    dict_cell_line['H58'] = {'genome' : 'hsa', 'tissue_source': 'kidney',   'cell_line':'A-704'}


    ## BaekLab - mRNA 6 Mouse  - mRNA-Seq
    dict_cell_line['H53'] = {'genome' : 'mmu', 'tissue_source': 'breast',     'cell_line':'MTV/TM-011'}
    dict_cell_line['H55'] = {'genome' : 'mmu', 'tissue_source': 'liver',      'cell_line':'Hepa-1c1c7'}
    dict_cell_line['H56'] = {'genome' : 'mmu', 'tissue_source': 'lung',       'cell_line':'LA-4'}
    dict_cell_line['H59'] = {'genome' : 'mmu', 'tissue_source': 'fibroblast', 'cell_line':'NOR-10'}
    dict_cell_line['H60'] = {'genome' : 'mmu', 'tissue_source': 'kidney',     'cell_line':'TCMK-1'}
    dict_cell_line['H61'] = {'genome' : 'mmu', 'tissue_source': 'kidney',     'cell_line':'Renca'}

    return dict_cell_line[sFileName]
#def END: assign_cell_line_info


def load_87_conserved_mirna_seedseq ():
    print('Loading 87 Conversed miRNA Seed Sequence')
    dict_sSeedSeq = {}
    InFile        = open('%s/87_conserved_mirna_full.txt' % sBASE_DIR, 'r')

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0           | 1       |
        # Column Descriptio: | miRNA(s0    | 7mer-m8 |
        # Column Example:    | hsa-let-7a; | GAGGTAG |

        list_sColumn = sReadLine.strip('\n').split('\t')

        # V-S Check:
        if len(list_sColumn) != nCONSERVED_MIRNA_COLUMN_SIZE:
            sys.exit('ERROR: Conserved MiRNA Column Size= %d\tContent= %s' % (len(list_sColumn), list_sColumn))

        list_sMirName   = list_sColumn[0].split(';')
        s7m8Seq         = reverse_complement(list_sColumn[1])
        s7A1Seq         = s7m8Seq[1:] + 'A'
        s8merSeq        = s7m8Seq + 'A'
        s6merSeq        = s7m8Seq[1:]

        for sMirName in list_sMirName:

            if sMirName not in dict_sSeedSeq:
                dict_sSeedSeq[sMirName] = {}

            dict_sSeedSeq[sMirName] = {'8mer':s8merSeq, '7m8':s7m8Seq, '7A1':s7A1Seq, '6mer':s6merSeq}
        #loop END: sMirName
    #loop END: sReadLine
    InFile.close()

    #V-S Check: Dictionary Size
    if len(dict_sSeedSeq) != nCONSERVED_MIRNAS:
        sys,exit('ERROR: dict_sSeedSeq Size= %d' % len(dict_sSeedSeq))

    print('Loading 87 Conversed miRNA Seed Sequence...DONE')

    return dict_sSeedSeq
#def END: load_87_conserved_mirna_seedseq


def output_text (sOutDir, sMirName, list_cSegs):

    OutFile = open('%s/%s.txt' % (sOutDir, sMirName), 'w')
    sHeader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('ChrID', 'Strand', 'miRNA', 'SeedSite',
                                                        'StartPos', 'EndPos', 'SeedCons', 'dEnergy',
                                                        'lAU', 'Seqment', 'DotBracket')
    OutFile.write(sHeader)

    for cSeg in list_cSegs:
        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
               % (cSeg.sChrID, cSeg.sStrand, cSeg.sMirName, cSeg.sSeedSeq, cSeg.nStartPos,
                  cSeg.nEndPos, cSeg.fSeedCons, cSeg.fMFE, cSeg.fLocalAU, cSeg.sSegmentSeq, cSeg.sDotBracket)
        OutFile.write(sOut)
    #loop END: sChrID, sMirName, sSeedSeq, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket

    OutFile.close()
#def END: output_text


def output_text_v2 (sOutDir, sMirName, list_sOutput):

    OutFile = open('%s/%s.txt' % (sOutDir, sMirName), 'w')
    sHeader = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('ChrID', 'miRNA', 'SeedSite', 'Start,End', 'dEnergy', 'lAU', 'Seqment', 'DotBracket')

    OutFile.write(sHeader)

    for sChrID, sMirName, sSeedSeq, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket in list_sOutput:
        sOut = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (sChrID, sMirName, sSeedSeq, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket)
        OutFile.write(sOut)
    #loop END: sChrID, sMirName, sSeedSeq, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket

    OutFile.close()
#def END: output_text


def parse_GFF3_files (sAnalysis):
    print('Parsing GFF3 File', time.ctime())

    # Step 1: Load GFF3 file
    sInDir      = dict_sMIRBASE_DIR['hsa'] if sAnalysis in list_sCancerTypes or \
                                              sAnalysis == 'test' else dict_sMIRBASE_DIR[sAnalysis]
    InFile      = open(sInDir, 'r')

    # Step 2: Read and parse GFF3 file
    list_cGFF3  = cMirBaseData_read_GFF3_file(InFile)

    # V-S Check: Empty List Check
    if not list_cGFF3:
        sys.exit('ERROR: cGFF3 List Size= %d' % len(list_cGFF3))

    print('Parsing GFF3 File...DONE', time.ctime())
    return list_cGFF3
#def END: parse_GFF3_files


def plot_SS_on_SVG_compiled (sOutDir, list_cSeg, sDirLabel, nTop, nWindowSize):
    print('Ploting SVG')

    sOutDir     = '%s/SVGPlot/%s' % (sOutDir, sDirLabel)
    os.makedirs(sOutDir, exist_ok=True)

    sBuffer     = ' ' * nWindowSize

    # Set window size
    nBase_width  = 16
    nBase_height = 550
    nSvg_width   = 20 + (105 * 12) + 300
    nSvg_height  = 650 * nTop

    nText_height      = 0
    nSS_height        = 10

    list_sOutput = []
    list_sOutput.append('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%d" height="%d">\n'
                        % (nSvg_width, nSvg_height))
    list_sOutput.append('<rect width="%.2f" height="%.2f" style="fill:#efefef;"/>\n'
                        % (nSvg_width, nSvg_height))
    list_sOutput.append('<g transform="translate(%.3f,%.3f)" font-family="Arial">\n'
                        % (20, 40))

    OutFile = open('%s/positionlist_100bps.txt' % sBASE_DIR, 'w')

    for cSeg in list_cSeg[:nTop]:

        sMatchedSeq, fConserved = cSeg.list_sCons[0]

        sHeader = '%s%s : %s : %s:%s : Start,End= %s : mfe= %0.2f : localAU= %0.4f : Conserved= %0.4f' \
                    % (cSeg.sChrID, cSeg.sStrand, cSeg.sMirName, cSeg.sSeedType, cSeg.sSeedSeq, cSeg.nStartPos,
                       cSeg.fMFE, cSeg.fLocalAU, fConserved)


        nSegmentLen             = len(cSeg.sSegmentSeq)
        sSeedSeq_buffered       = sBuffer + cSeg.sSeedSeq + sBuffer

        sSeedSeq_buffered, nWindowStart, nWindowEnd = mark_additional_paired_bps (cSeg, sSeedSeq_buffered, nWindowSize, 0)

        OutFile.write('%s\t%s:%s-%s\n' % (cSeg.sMirName, cSeg.sChrID,nWindowStart, nWindowEnd))

        #print(cSeg.sChrID, cSeg.sMirName, cSeg.list_nPos)
        #print(cSeg.sSegmentSeq)
        #print(len(cSeg.sSegmentSeq))
        #print(cSeg.sDotBracket)
        #print(len(cSeg.sDotBracket))
        #print(sSeedSeq_buffered)
        #print(len(sSeedSeq_buffered))
        #print()

        list_sOutput.append('<text x="20" y="%s" font-size="20" style="font-weight:bold;">MutPanel: %s</text>\n'
                            % (nText_height,sHeader))

        sOutput, nSS_max_height = JM_rnastructure(cSeg.sSegmentSeq, cSeg.sDotBracket,
                                                  sSeedSeq_buffered, nSegmentLen,(nSvg_width/2) - 150,
                                                  nSS_height)

        if nSS_max_height > nBase_height:
            nHeight_increase = nSS_max_height - nBase_height
        else:
            nHeight_increase = 0

        nText_height            = nText_height +  nBase_height + nHeight_increase
        nSS_height              = nSS_height   +  nBase_height + nHeight_increase

        list_sOutput.append(sOutput)  ## rna structure - normal sample
    #loop END: cSeg

    OutFile.close()

    list_sOutput.append('</g>\n')
    list_sOutput.append('</svg>\n')

    sOutput = ''.join(list_sOutput)

    sOutFileDir  = '%s/%s.%s.%s.svg' % (sOutDir, cSeg.sChrID+cSeg.sStrand, cSeg.sMirName,cSeg.nStartPos)
    OutFile2     = open(sOutFileDir, 'w')
    OutFile2.write(sOutput)
    OutFile2.close()
    ## [3] convert .svg into .png file
    svg2png(sOutFileDir)

    print('Ploting SVG...DONE')
#def END: plot_SS_on_SVG_compiled


def plot_SS_on_SVG_individual (sOutDir, list_cSeg, sDirLabel, nWindowSize, sTargetSeedType, sSortBy):
    print('Ploting SVG')

    sOutDir     = '%s/SVGPlot/%s' % (sOutDir, sDirLabel)
    os.makedirs(sOutDir, exist_ok=True)

    OutFile = open('%s/%s_positionlist_%sbps%scons_%s.txt' % (sBASE_DIR, sTargetSeedType, nWindowSize, nSECOND_CONS_FILTER, sSortBy), 'w')

    nCnt = 0


    for cSeg in list_cSeg[:100]:

        print(cSeg.sMirName)

        nSS_height               = 80   # x coordinate of secondary structure

        nSegmentLen              = len(cSeg.sSegmentSeq)

        sSeedSeq_buffered        = ' ' * cSeg.nLowerWin + cSeg.sSeedSeq + ' ' * cSeg.nUpperWin

        #V-S Check: Sequence Size
        if nSegmentLen != len(sSeedSeq_buffered):
            sys.exit('ERROR: Invalid Size\tSegment= %d\tBuffered= %d' % (nSegmentLen, len(sSeedSeq_buffered)))

        sSeedSeq_buffered        = mark_additional_paired_bps (cSeg, sSeedSeq_buffered, 0)

        list_sOutput = get_SVG_outputlist (cSeg)

        OutFile.write('%s\t%0.3f\t%0.3f\t%0.3f\t%s:%s-%s\n'
                      % (cSeg.sMirName, cSeg.fMFE, cSeg.fLocalAU,
                         cSeg.fConsScore, cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd))


        sOutput, nSS_max_height = JM_rnastructure(cSeg.sSegmentSeq, cSeg.sDotBracket,
                                                  sSeedSeq_buffered, nSegmentLen, 150,
                                                  nSS_height)


        list_sOutput.append(sOutput)  ## rna structure - normal sample


        list_sOutput.append('</g>\n')
        list_sOutput.append('</svg>\n')

        sOutput = ''.join(list_sOutput)

        sOutFileDir  = '%s/%s_%s.%s.%s.svg' % (sOutDir, nCnt+1, cSeg.sChrID+cSeg.sStrand, cSeg.sMirName,cSeg.nStartPos)
        OutFile2     = open(sOutFileDir, 'w')
        OutFile2.write(sOutput)
        OutFile2.close()
        ## [3] convert .svg into .png file
        svg2png(sOutFileDir)

        nCnt += 1
    #loop END: cSeg
    OutFile.close()
    print('Ploting SVG...DONE')

    return list_cSeg
#def END: plot_SS_on_SVG_individual


def plot_SS_on_SVG_extended (sOutDir, list_cSeg, sDirLabel, nWindowSize, sTargetSeedType, sSortBy):
    print('Ploting SVG')

    sOutDir         = '%s/SVGPlot/%s' % (sOutDir, sDirLabel)
    os.makedirs(sOutDir, exist_ok=True)

    nCnt            = 0

    OutFile = open('%s/%s_positionlist_%sbps%scons_%s_ext.txt' % (sBASE_DIR, sTargetSeedType, nWindowSize, nSECOND_CONS_FILTER, sSortBy), 'w')

    for cSeg in list_cSeg[:100]:
        nSvg_width              = 20 + (105 * 12) + 300
        nSS_height              = 80

        cSeg, bOpenSeedFlag     = extend_segment_info (cSeg, nWindowSize, nEXTEND_LEN, sTargetSeedType)


        if bOpenSeedFlag == False: continue

        # original
        nSegmentLen              = len(cSeg.sSegmentSeq)
        sSeedSeq_buffered        = ' ' * cSeg.nLowerWin + cSeg.sSeedSeq + ' ' * cSeg.nUpperWin

        #V-S Check: Sequence Size
        if nSegmentLen != len(sSeedSeq_buffered):
            sys.exit('ERROR: Invalid Size\tSegment= %d\tBuffered= %d' % (nSegmentLen, len(sSeedSeq_buffered)))


        list_sOutput = get_SVG_outputlist (cSeg)


        sSeedSeq_buffered = mark_additional_paired_bps (cSeg, sSeedSeq_buffered, 0)

        OutFile.write('%s\t%0.3f\t%0.3f\t%0.3f\t%s:%s-%s\n'
                      % (cSeg.sMirName, cSeg.fMFE, cSeg.fLocalAU,
                         cSeg.fConsScore, cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd))

        sOutput, nSS_max_height = JM_rnastructure(cSeg.sSegmentSeq, cSeg.sDotBracket,
                                                  sSeedSeq_buffered, nSegmentLen, 0,
                                                  nSS_height)

        list_sOutput.append(sOutput)


        # Extended
        nSegmentLen             = len(cSeg.sSegmentSeq_ext)
        sSeedSeq_buffered       = ' ' * (cSeg.nLowerWin + nEXTEND_LEN) + cSeg.sSeedSeq + ' ' * (cSeg.nUpperWin + nEXTEND_LEN)

        #V-S Check: Sequence Size
        if nSegmentLen != len(sSeedSeq_buffered):
            sys.exit('ERROR: Invalid Size\tSegment= %d\tBuffered= %d' % (nSegmentLen, len(sSeedSeq_buffered)))


        sSeedSeq_buffered_ext   = mark_additional_paired_bps (cSeg, sSeedSeq_buffered, nEXTEND_LEN)


        sOutput, nSS_max_height = JM_rnastructure(cSeg.sSegmentSeq_ext, cSeg.sDotBracket_ext,
                                                  sSeedSeq_buffered_ext, nSegmentLen,((nSvg_width)/2)-50,
                                                  nSS_height)

        list_sOutput.append(sOutput)  ## rna structure - normal sample

        list_sOutput.append('</g>\n')
        list_sOutput.append('</svg>\n')

        sOutput = ''.join(list_sOutput)

        sOutFileDir  = '%s/%s_%s.%s.%s.svg' % (sOutDir, nCnt+1, cSeg.sChrID+cSeg.sStrand, cSeg.sMirName,cSeg.nStartPos)
        OutFile2      = open(sOutFileDir, 'w')
        OutFile2.write(sOutput)
        OutFile2.close()
        ## [3] convert .svg into .png file
        svg2png(sOutFileDir)

        nCnt += 1
    #loop END: cSeg
    OutFile.close()

    print('Ploting SVG...DONE')
#def END: plot_SS_on_SVG_v2


def get_SVG_outputlist (cSeg):
    # Set window size
    nSvg_width   = 20 + (105 * 12) + 300
    nSvg_height  = 1000

    nText_height      = 0

    list_sOutput = []
    list_sOutput.append('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%d" height="%d" xml:space="preserve">\n'
                        % (nSvg_width, nSvg_height))
    list_sOutput.append('<rect width="%.2f" height="%.2f" style="fill:#efefef;"/>\n'
                        % (nSvg_width, nSvg_height))
    list_sOutput.append('<g transform="translate(%.3f,%.3f)" font-family="Arial">\n'
                        % (20,40))

    fConserved = cSeg.fConsScore


    sHeader = '%s%s   :   %s   :   Site= %s:%s   :   InterDist= %d  :   AddPair= %smer:%s   :   Start,End= %s,%s' % \
              (cSeg.sChrID, cSeg.sStrand, cSeg.sMirName, cSeg.sSeedType,
               cSeg.sSeedSeq, cSeg.nInterDist, len(cSeg.sAddPair), cSeg.sAddPair, cSeg.nStartPos, cSeg.nEndPos)

    list_sOutput.append('<text x="20" y="%s" font-size="20" style="font-weight:bold;">MutPanel: %s</text>\n'
                             % (nText_height,sHeader))


    sHeader1 = 'mfe              = %0.2f'   % cSeg.fMFE
    sHeader2 = 'localAU        = %0.4f'     % cSeg.fLocalAU
    sHeader3 = 'conservation= %0.4f'        % fConserved

    list_sOutput.append('<text x="20" y="%s" font-size="20" style="font-weight:bold;">                %s</text>\n'
                        % (nText_height + 20,sHeader1))

    list_sOutput.append('<text x="20" y="%s" font-size="20" style="font-weight:bold;">                %s</text>\n'
                        % (nText_height + 40,sHeader2))

    list_sOutput.append('<text x="20" y="%s" font-size="20" style="font-weight:bold;">                %s</text>\n'
                        % (nText_height + 60,sHeader3))
    return list_sOutput
#def END: get_SVG_outputlist


def mark_additional_paired_bps (cSeg, sSeedSeq_buffered, nEXTEND_LEN):

    list_sSeq = list(sSeedSeq_buffered)

    for sOpenInfo in cSeg.list_sOpenSeq:

        nStartPos = int(sOpenInfo.split(',')[0]) + nEXTEND_LEN
        nEndPos   = int(sOpenInfo.split(',')[1]) + nEXTEND_LEN

        sOpenSeq    = sOpenInfo.split(',')[2]

        for sReIndex in re.finditer(reverse_complement(cSeg.list_sAdd_Cons[0][0]), sOpenSeq):

            nReIndexStart = sReIndex.start()
            nReIndexEnd   = sReIndex.end()

            nIndexStart = nStartPos   + nReIndexStart
            nIndexEnd   = nStartPos   + nReIndexEnd

            ## For Output Check on UCSC Genome Browser
            if cSeg.sStrand == '+':

                nWindowStart = int(cSeg.nStartPos) + nIndexStart  + 1
                nWindowEnd   = int(cSeg.nStartPos) + (cSeg.nLowerWin + len(cSeg.sSeedSeq))

            else:
                nWindowStart = int(cSeg.nEndPos) - nIndexStart  + 1
                nWindowEnd   = int(cSeg.nEndPos) - (cSeg.nLowerWin + len(cSeg.sSeedSeq))

            #if END: cSeg.sStrand

            sTargetSeq        = cSeg.sSegmentSeq[nIndexStart:cSeg.nSeedEnd] #Seed + InterDist + AddPair
            cSeg.sAddPair     = cSeg.sSegmentSeq[nIndexStart:nIndexEnd]
            cSeg.nInterDist   = len(sTargetSeq) - len(cSeg.sSeedSeq) - len(cSeg.sAddPair)
            cSeg.nWindowStart = nWindowStart
            cSeg.nWindowEnd   = nWindowEnd

            for i in range(nIndexStart, nIndexEnd):
                list_sSeq[i] = sOpenSeq[i - nStartPos]
            #loop END: i
        #loop END: sReIndex
    #loop END: sOpenInfo

    sNewSeq = ''.join(list_sSeq)

    return sNewSeq
#def END: mark_additional_paired_bps


def JM_rnastructure(sSegmentSeq, sDotBracket, sBufferedSeedSeq, sSegmentLen, gx, gy):

    list_sOutput = []
    list_sOutput.append('<g transform="translate(%.3f,%.3f)">\n' % (gx, gy))

    ss_pairidx = pairidx(sDotBracket, 2)
    ss_coor    = rnaplot(sSegmentSeq, sDotBracket) # Coordintes
    plp = []

    for i in range(sSegmentLen):
        plx, ply = ss_coor[i]
        plp.append('%.2f,%.2f'%(plx, ply))
    #loop END: i
    list_sOutput.append('<polyline points="%s" style="fill:none;stroke:black;stroke-width:1;"/>\n'
                        %   (' '.join(plp)))

    for i, j in ss_pairidx:
        lx1 = ss_coor[i][0]
        ly1 = ss_coor[i][1]
        lx2 = ss_coor[j][0]
        ly2 = ss_coor[j][1]
        list_sOutput.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" style="stroke:black;"/>\n'
                            %   (lx1, ly1, lx2, ly2))
    #loop END: i, j

    nFontSize = 10
    nRadius   = 5


    for i in range(sSegmentLen):
        nCircle_X, nCircle_Y  = ss_coor[i] # Circle X and Y coordinate

        sCircleFill   = 'white' if sBufferedSeedSeq[i] == ' ' else 'black'  # Circle Fill Color
        cstroke       = 'black'

        list_sOutput.append('<circle cx="%.2f" cy="%.2f" r="%.2f" style="fill:%s;stroke:%s;"/>\n'
                            %   (nCircle_X, nCircle_Y, nRadius, sCircleFill, cstroke))

        nText_X = nCircle_X
        nText_Y = nCircle_Y + nFontSize * 0.35

        sTextColor = 'white' if sSegmentSeq[i] == sBufferedSeedSeq[i] else 'black'
        list_sOutput.append('<text x="%.2f" y="%.2f" font-size="%.2f" text-anchor="middle" style="fill:%s;">%s</text>\n'
                            %   (nText_X, nText_Y, nFontSize, sTextColor, sSegmentSeq[i]))
    #loop END: i

    list_sOutput.append('</g>\n')

    return ''.join(list_sOutput), max([nCircle_Y for nCircle_X, nCircle_Y in ss_coor])
#def END: JM_rnastructure

## endregion


## region Classes

## region class cFasta
re_nonchr = re.compile('[^a-zA-Z]')
class cFasta:
    def __init__(self, sRefFile):

        #V-S Check: File Existence
        if not os.path.isfile(sRefFile):
           sys.exit('(): File does not exist')

        self.InFile     = open(sRefFile, 'r')
        self.sChrIDList = []
        self.nChromLen  = []
        self.nSeekPos   = []
        self.nLen1      = []
        self.nLen2      = []

        #V-S Check: File Existence
        if not os.path.isfile('%s.fai'%sRefFile):
           sys.exit('.fai file does not exist')

        InFile = open('%s.fai' % sRefFile, 'r')
        for sLine in InFile:
           list_sColumn = sLine.strip('\n').split() # Goes backwards, -1 skips the new line character

           self.sChrIDList.append  (list_sColumn[0])
           self.nChromLen.append   (int(list_sColumn[1]))
           self.nSeekPos.append    (int(list_sColumn[2]))
           self.nLen1.append       (int(list_sColumn[3]))
           self.nLen2.append       (int(list_sColumn[4]))
        #loop END: sLINE
        InFile.close()
        self.sType = []
    #def END: __init_

    def fetch(self, sChrom, nFrom = None, nTo = None, sStrand = '+'):
        assert sChrom in self.sChrIDList, sChrom
        nChrom = self.sChrIDList.index(sChrom)

        if nFrom == None: nFrom = 0
        if nTo   == None: nTo = self.nChromLen[nChrom]

        assert(0 <= nFrom) and(nFrom < nTo) and(nTo <= self.nChromLen[nChrom])

        nBlank = self.nLen2[nChrom] - self.nLen1[nChrom]

        nFrom  = int(nFrom +(nFrom / self.nLen1[nChrom]) * nBlank) # Start Fetch Position

        nTo    = int(nTo   +(nTo   / self.nLen1[nChrom]) * nBlank) # End Fetch Position

        self.InFile.seek(self.nSeekPos[nChrom] + nFrom)            # Get Sequence

        sFetchedSeq = re.sub(re_nonchr, '', self.InFile.read(nTo - nFrom)).upper()

        if   sStrand == '+':
            return sFetchedSeq

        elif sStrand == '-':
            return reverse_complement(sFetchedSeq)

        else:
            sys.exit('Error: moudle seq: invalid strand')
        #if END: sStrand
    #def END: fetch
#class END: Fasta
## endregion

## region class cSegment
class cSegment:
    def __init__(self):
        self.sChrID      = 'NULL'
        self.sStrand     = 'NULL'
        self.sMirName    = 'NULL'
        self.sSeedType   = 'NULL'
        self.nSeedStart  = 0
        self.nSeedEnd    = 0
        self.sSeedSeq    = 'NULL'
        self.fSeedCons   = 0.0
        self.nStartPos   = 0
        self.nEndPos     = 0
        self.fMFE        = 0.0
        self.fLocalAU    = 0.0
        self.sSegmentSeq = 'NULL'
        self.sDotBracket = 'NULL'
    #def END: __init__

def cSegment_assign_segment_data (cSeg, sChrID, sMirName, sSeedType, nSeedStart, nSeedEnd,  sSeedSeq, sPosition, fMFE, fLocalAU, fSeedCons, sSegmentSeq, sDotBracket):

    cSeg.sChrID      = sChrID[:-1]
    cSeg.sStrand     = sChrID[-1]
    cSeg.sMirName    = sMirName
    cSeg.sSeedType   = sSeedType
    cSeg.nSeedStart  = nSeedStart
    cSeg.nSeedEnd    = nSeedEnd
    cSeg.sSeedSeq    = sSeedSeq
    cSeg.nStartPos   = int(sPosition.split(',')[0])
    cSeg.nEndPos     = int(sPosition.split(',')[1])
    cSeg.fMFE        = fMFE
    cSeg.fLocalAU    = fLocalAU
    cSeg.fSeedCons   = fSeedCons
    cSeg.sSegmentSeq = sSegmentSeq
    cSeg.sDotBracket = sDotBracket

    return cSeg
#def END: cSegment_assign_segment_data


def cSegment_assign_segment_data_v2 (cSeq, list_sSegmentData):
    # sSegmentData Index Format:
    # Index:    0       | 1           | 2              | 3                       | 4       | 5                               |
    # Data:     sChrID  | sMirName    | sSeedInfo      | nPos                    | fMFE | fLocalAu                        |
    # Example:  chr1+   |  miR-99a-5p | 8mer:TACGGGTA  | ['32079366,32079474']   | -25.2   |  0.53177848                     |
    # Index:    6                                                                                                            |
    # Data:     sSegmentSeq                                                                                                  |
    # Example:  CTCCGTCTCTGCAGCCTACCCAGCCTCCGAGTTTGTTTGTTGAGTGAGTATACGGGTAAGTGGACAGACAGAAGGCTGCCAGCTCCAGCCCCTTCTTCCTGGCTCCTG |
    # Index:    7                                                                                                            |
    # Data:     sDotBracket                                                                                                  |
    # Example:  ........(((.(((.....((((((.....((((((((((.(.(.............).).))))))))))))))))...))).)))..((........))...... |
    # Index:    8                   | 9
    # Data:     list_sOpen          | list_sCons
    # Example:  ['16,24,AGCACTTA']  | [['AAGTGCT', 59.761001110076904]]]




    #V-S Check: Index Number
    if len(list_sSegmentData) != 10:
        sys.exit('ERROR: list_sSegmentData Size= %d \n Contents: \n %s'
                 % (len(list_sSegmentData), list_sSegmentData))

    cSeq.sChrID      = list_sSegmentData[0][:-1]
    cSeq.sStrand     = list_sSegmentData[0][-1]
    cSeq.sMirName    = list_sSegmentData[1]
    cSeq.sSeedType   = list_sSegmentData[2].split(':')[0]
    cSeq.sSeedSeq    = list_sSegmentData[2].split(':')[1]
    cSeq.list_nPos   = list_sSegmentData[3]
    cSeq.fMFE     = list_sSegmentData[4]
    cSeq.fLocalAU    = list_sSegmentData[5]
    cSeq.sSegmentSeq = list_sSegmentData[6]
    cSeq.sDotBracket = list_sSegmentData[7]
    cSeq.list_sOpen  = list_sSegmentData[8]
    cSeq.list_sCons  = list_sSegmentData[9]
#def END: cSegment_assign_segment_data
#class END: Segment



## endregion

## region class cMicroRNAData

class cMicroRNAData:
    def __init__(self):
        self.sGenus         = 'NULL'
        self.sSpecies       = 'NULL'
        self.sMiRNASeq      = 'NULL'
        self.sMiRNASym      = 'NULL'
        self.sMiRNASeedSeq  = 'NULL'
    #def END: __init__


    def parse_miRNA_sym(self, sReadLine):
        sMirTitle, sMirID, sGenus, sSpecies, sMiRNASym = sReadLine.split(' ')

        self.sGenus         = sGenus
        self.sSpecies       = sSpecies
        self.sMiRNASym      = sMiRNASym.replace('\n', '')
    #def END: parse_miRNA_sym


    def parse_miRNA_seq(self, sReadLine):
        self.sMiRNASeq      = sReadLine.replace('\n', '')
    #def END: parse_miRNA_seq


    def determine_seed_seq(self):
        sBaseComplementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}

        sBasesList          = list(self.sMiRNASeq[0:8]) # Turns the sequence in to a list

        sBasesList          = [sBaseComplementDict[sBase] for sBase in sBasesList]

        self.sMiRNASeedSeq  = ''.join(sBasesList)[::-1]
        # For each base, the dictionary value for that base key is assigned
    #def END: determine_seed_seq
#class END: cMicroRNAData

## endregion

## region class cMirBaseData
class cMirBaseData:
    def __init__(self):
        self.sGenome    = 'NULL'
        self.sMirName   = 'NULL'
        self.sChrID     = 'NULL'
        self.sStrand    = 'NULL'
        self.sFeature   = 'NULL'
        self.nStartPos  = 0
        self.nEndPos    = 0
    #def END: __init__
#class END: cMirBaseData


def cMirBaseData_read_GFF3_file (InFile):
    # GFF3 File Format
    # Column Number:        0             | 1       | 2                         | 3            | 4            | 5             |
    # Column Description:   contig        | SPACE   | feature                   | start pos    | end pos      | SPACE         |
    # Column Example:       chr1          | .       | miRNA_primary_transcript  | 30366        | 30503        | .             |
    #
    # Column Number:        6             | 7       | 8                                                                       |
    # Column Description:   strand        | SPACE   | attributes                                                              |
    # Column Example:       +             | .       | ID=MI0006363_1;accession_number=MI0006363;Name=hsa-mir-1302-2           |
    list_cGFF3 = []
    for sReadLine in InFile:
        cGFF3 = cMirBaseData()

        if sReadLine.startswith('#'): continue
        else:
            list_sColumn = sReadLine.strip('\n').split('\t')

            # V-S Check: Column List Size
            if len(list_sColumn) != nGFF3_COLUMN_SIZE:
                sys.exit('ERROR: GFF3 Column List Size= %d' % len(list_sColumn))

            cGFF3.sChrID     = list_sColumn[0] if 'chr' in list_sColumn[0] else 'chr%s' % list_sColumn[0]
            cGFF3.sFeature   = list_sColumn[2]
            cGFF3.nStartPos  = int(list_sColumn[3])  # 1-based closed
            cGFF3.nEndPos    = int(list_sColumn[4])  # 1-based closed
            cGFF3.sStrand    = list_sColumn[6]

            sAttribute       = list_sColumn[8].split(';')
            cGFF3.sMirName   = sAttribute[2].split('=')[1].lower()
            cGFF3.sGenome    = cGFF3.sMirName[:3]

            list_cGFF3.append(cGFF3)
        #if END: sReadLine.startswith

    #loop END: sReadLine

    # V-S Check: Empty List Check
    if not list_cGFF3:
        sys.exit('ERROR: cGFF3 List Size= %d' % len(list_cGFF3))

    return list_cGFF3
#class END: cMirBaseData_read_GFF3_file
## endregion

## endregion


## region Segmentation

def segmentation (sOutDir, sGenome, sSeedType, sChrID, sStrand, nWindowSize):
    print('Segmentation Processing %s' % time.ctime())

    cGenome         = cFasta(dict_sGENOME_DIR[sGenome])
    dict_sSeedSeq   = load_87_conserved_mirna_seedseq()

    list_sMirName   = sorted(list(dict_sSeedSeq.keys()))

    sChrSeq         = cGenome.fetch(sChrID).upper()

    for sMirName in list_sMirName:

        sChrDir = '%s/%s%s' % (sOutDir, sChrID, sStrand)
        os.makedirs (sChrDir, exist_ok=True)

        # File Check:
        if os.path.isfile('%s/%s.segments' % (sChrDir, sMirName)):
            print('%s Processed' % sMirName)
            continue

        print('Processing ---%s---  *** %s%s ***' % (sMirName, sChrID, sStrand))

        if sStrand == '+': sSeedSeq = dict_sSeedSeq[sMirName][sSeedType]
        else:              sSeedSeq = reverse_complement(dict_sSeedSeq[sMirName][sSeedType])

        dict_sSegment = get_segment (sSeedSeq, sChrSeq, sStrand, int(nWindowSize))

        pickle.dump(dict_sSegment, open('%s/%s.segments' % (sChrDir, sMirName), 'wb'))
    #loop END: sMirName

    print('Segmentation Processing...DONE %s' % time.ctime())
#def END: segmentation


def get_segment (sSeedSeq, sChrSeq, sStrand, nWindowSize):
    print('Getting Segments')

    dict_sSegment = {}

    nCntwN        = 0

    for sReIndex in re.finditer(sSeedSeq, sChrSeq):
        nIndexStart = sReIndex.start()
        nIndexEnd   = sReIndex.end()

        nStartPos   = nIndexStart - nWindowSize
        nEndPos     = nIndexEnd   + nWindowSize

        if nStartPos < 0:
            print(nStartPos, nEndPos)
            continue

        if sStrand == '+':  sSegmentSeq = sChrSeq[nStartPos:nEndPos]
        else:               sSegmentSeq = reverse_complement(sChrSeq[nStartPos:nEndPos])

        if 'N' in sSegmentSeq: nCntwN += 1
        else:

            sPositions = '%s,%s' % (nStartPos, nEndPos)

            if sPositions not in dict_sSegment:
                dict_sSegment[sPositions] = ''
            dict_sSegment[sPositions] = sSegmentSeq
        #if END: 'N'
    #loop END: sReIndex

    print('Getting Segments...DONE')
    print('Segments', len(dict_sSegment))
    print('WIth N', nCntwN)

    #V-S Check: Empty Dictionary
    #if not dict_sSegment:
        #sys.exit('ERROR: dict_sSegment Size= %d' % len(dict_sSegment))

    return dict_sSegment
#def END: get_segment

## endregion


## region Filter: Seed Conservation

def filter_by_seed_conservation (sInDir, sOutDir, sGenome, sChrID, sSeedType, nWindowSize):

    print('Filtering By SeedSite Conversation Filter= %d %s' % (nFIRST_CONS_FILTER, time.ctime()))

    sGenomeDir      = 'hg19' if sGenome == 'hsa' else 'mm9'
    sVertDir        = 'vertebrate_100way' if sGenome == 'hsa' else 'vertebrate'

    nWindowSize     = int(nWindowSize)
    dict_sSeedSeq   = load_87_conserved_mirna_seedseq ()
    nSeedLength     = dict_nSEED_LENGTH[sSeedType]

    list_sMirName   = sorted(list(dict_sSeedSeq.keys()))

    for sMirName in list_sMirName:

        dict_sSegement_filtered = {}

        # File Check:
        if os.path.isfile('%s/%s.filt.segments' % (sOutDir, sMirName)):
            print('%s Processed' % sMirName)
            continue

        dict_sSegment = pickle.load(open('%s/%s/%s.segments' % (sInDir, sChrID, sMirName), 'rb'))

        for sPosition in dict_sSegment:

            nSeedStartPos   = int(sPosition.split(',')[0]) + nWindowSize
            nSeedEndPos     = int(sPosition.split(',')[1]) - nWindowSize
            sSegmentSeq     = dict_sSegment[sPosition]

            list_fCons      = fetch_cbc_cons(sGenomeDir, sVertDir, sChrID[:-1], nSeedStartPos, nSeedEndPos)

            #Filter by Seed Site Conservation Score
            if sum(list_fCons) < nFIRST_CONS_FILTER: continue

            if sPosition not in dict_sSegement_filtered:
                dict_sSegement_filtered[sPosition] = []

            dict_sSegement_filtered[sPosition] = [sSegmentSeq, sum(list_fCons)]
        #loop END: sPosition

        print('%s Before= %d \t After= %d' % (sMirName, len(dict_sSegment), len(dict_sSegement_filtered)))

        pickle.dump(dict_sSegement_filtered, open('%s/%s.segments' % (sOutDir, sMirName), 'wb'))
    #loop END: sMirName
    print('Filtering By SeedSite Conversation Filter= %d...DONE %s' % (nFIRST_CONS_FILTER, time.ctime()))
#def END: filter_by_seed_conservation

## endregion


## region RNAfold and locaAU

def RNAfold_and_localAU (sInDir, sOutDir, sGenome, sChrID, sSeedType, nWindowSize):
    print('Processing RNAfold and LocalAU %s' % time.ctime())

    dict_sSeedSeq   = load_87_conserved_mirna_seedseq()
    nSeedLength     = dict_nSEED_LENGTH[sSeedType]

    list_sMirName   = sorted(list(dict_sSeedSeq.keys()))

    for sMirName in list_sMirName:
        #sMirName = 'miR-214-3p'

        #V-S Check: File Path
        if os.path.isfile('%s/%s.list' % (sOutDir, sMirName)):
            print('%s Processed' % sMirName)
            continue

        print('Processing %s %s' % (sMirName, time.ctime()))
        list_cSegs = []

        # Dictionary Key = Start,Pos    Value = [sSegmentSeq, SeedCons]
        dict_sSegment  = pickle.load(open('%s/%s/%s.segments' % (sInDir, sChrID, sMirName), 'rb'))

        for sPosition in dict_sSegment:

            cSeg = cSegment()

            sSegmentSeq = dict_sSegment[sPosition][0]
            fSeedCons   = dict_sSegment[sPosition][1]

            if '-' in cSeg.sSegmentSeq: continue

            # Get DotBracket Notation and Min. Free Energy : Code Credit - Sukjun Kim
            sDotBracket, fMFE = get_dotbracket_and_MFE (sSegmentSeq, 0.0)

            nSeedStartPos        = int(nWindowSize)
            nSeedEndPos          = int(nWindowSize) + nSeedLength

            # Get LocalAU Score : Code Credit - Doyeon Kim
            fLocalAU        = get_LocalAU_score (nSeedStartPos, nSeedEndPos, sSegmentSeq, nSeedLength, int(nWindowSize))

            cSeg = cSegment_assign_segment_data(cSeg, sChrID, sMirName, sSeedType, nSeedStartPos, nSeedEndPos,
                                                dict_sSeedSeq[sMirName][sSeedType],
                                                sPosition, fMFE, fLocalAU, fSeedCons, sSegmentSeq,
                                                sDotBracket )
            list_cSegs.append(cSeg)
        #loop END: sSegmentSeq

        list_cSegs = sorted(list_cSegs, key=lambda c: c.fSeedCons, reverse=True)

        # Text Output
        output_text (sOutDir, sMirName, list_cSegs)

        # Pickle Dump
        pickle.dump(list_cSegs, open('%s/%s.list' % (sOutDir, sMirName), 'wb'))

        print('Processing %s...DONE %s' % (sMirName, time.ctime()))
    #loop END: sMirName
    print('Processing RNAfold and LocalAU...DONE %s' % time.ctime())
#def END: RNAfold


def RNAfold_and_localAU_v2 (sInDir, sOutDir, sGenome, sChrID, sSeedType, nWindowSize):
    print('Processing RNAfold and LocalAU %s' % time.ctime())

    dict_sSeedSeq   = load_87_conserved_mirna_seedseq ()
    nSeedLength     = dict_nSEED_LENGTH[sSeedType]

    list_sMirName   = sorted(list(dict_sSeedSeq.keys()))

    for sMirName in list_sMirName:
        #sMirName = 'miR-214-3p'

        #V-S Check: File Path
        if os.path.isfile('%s/%s.list' % (sOutDir, sMirName)):
            print('%s Processed' % sMirName)
            continue

        print('Processing %s %s' % (sMirName, time.ctime()))
        list_cSegs = []

        dict_sSegment  = pickle.load(open('%s/%s/%s.segments' % (sInDir, sChrID, sMirName), 'rb'))

        for sPosition in dict_sSegment:

            sSegmentSeq = dict_sSegment[sPosition]

            if '-' in sSegmentSeq: continue

            # Get DotBracket Notation and Min. Free Energy : Code Credit - Sukjun Kim
            sDotBracket, fMFE = get_dotbracket_and_MFE (sSegmentSeq, 0.0)

            nSeedStartPos        = int(nWindowSize)
            nSeedEndPos          = int(nWindowSize) + nSeedLength

            # Get LocalAU Score : Code Credit - Doyeon Kim
            fLocalAUScore        = get_LocalAU_score (nSeedStartPos, nSeedEndPos, sSegmentSeq, nSeedLength, int(nWindowSize))

            list_cSegs.append([sChrID, sMirName, '%s:%s' % (sSeedType, dict_sSeedSeq[sMirName][sSeedType]),
                                 sPosition, fMFE, fLocalAUScore, sSegmentSeq, sDotBracket])
        #loop END: sSegmentSeq

        list_cSegs = sorted(list_cSegs, key=lambda e:(e[4], -e[5]))

        # Text Output
        output_text (sOutDir, sMirName, list_cSegs)

        # Pickle Dump
        pickle.dump(list_cSegs, open('%s/%s.list' % (sOutDir, sMirName), 'wb'))

        print('Processing %s...DONE %s' % (sMirName, time.ctime()))
    #loop END: sMirName
    print('Processing RNAfold and LocalAU...DONE %s' % time.ctime())
#def END: RNAfold


def get_dotbracket_and_MFE (sSegmentSeq,  fMFE):

    sScript      = 'echo "%s" | %s -s -e %d' % (sSegmentSeq, sRNAfold_DIR, fMFE)
    list_sStdOut = subprocess.Popen(sScript, stdout=subprocess.PIPE, shell=True).stdout

    list_sDotBr  = []

    for i, sReadLine in enumerate(list_sStdOut):
        sReadLine = str(sReadLine, 'UTF-8').strip('\n')

        if i >= 1:
            list_sColumn = sReadLine.split()
            sDotBracket  = list_sColumn[0]

            fMFE      = float(list_sColumn[1])
            list_sDotBr.append([sDotBracket, fMFE])
    #loop END: i, sReadLine

    return list_sDotBr[0]
#def END: get_dotbracket_and_MFE


def get_LocalAU_score(nStartPos, nEndPos, sSeq, nSeedLength, nWindowSize):

    list_nPosWeight = get_weighted_pos_list(nSeedLength, nWindowSize)

    if nStartPos < nWindowSize:
        nSeqOffset = nWindowSize - nStartPos
    else:
        nSeqOffset = 0

    list_nPosWeight = list_nPosWeight[nSeqOffset:]

    nSeqFrom        = nStartPos - nWindowSize + nSeqOffset
    nSeqTo          = nSeqFrom + len(list_nPosWeight)

    sSubSeq         = sSeq[nSeqFrom:nSeqTo]

    fAUScore        = 0.0
    fMaxValue       = 0.0

    for nSeqPos in range(0, len(sSubSeq)):

        if sSubSeq[nSeqPos] == 'A' or sSubSeq[nSeqPos] == 'U' or sSubSeq[nSeqPos] == 'T' :
            fAUScore += list_nPosWeight[nSeqPos]

        fMaxValue += list_nPosWeight[nSeqPos]
    #loop END: nSeqPos

    fAUScore /= fMaxValue
    return fAUScore
#def END: get_LocalAU


def get_weighted_pos_list (nSeedLength, nWindowSize):

    list_nPosWeight = []
    nStart          = 2

    for n in range(0, nWindowSize) :
        list_nPosWeight.append( float(1.0)/float(nStart+n)  )
    list_nPosWeight.reverse()

    for n in range(0, nSeedLength) :
        list_nPosWeight.append( 0 )

    for n in range(0, nWindowSize) :
        list_nPosWeight.append( float(1.0)/float(nStart+n) )

    return list_nPosWeight
#def END: get_localAU

## endregion


## region Screen Dot-Bracket - Part 1: Open Seed Site

def filter_dot_bracket(sInDir, sOutDir, sGenome,  sChrID, nWindowSize):
    print('Filtering Dot-Bracket %s' % time.ctime())

    dict_sSeedSeq = load_87_conserved_mirna_seedseq()
    list_sMirName = sorted(list(dict_sSeedSeq.keys()))
    list_cGFF3    = parse_GFF3_files(sGenome)

    for sMirName in list_sMirName:
        print('Filtering %s %s' % (sMirName, sChrID))

        if sMirName == 'miR-6766-3p': continue
        if sMirName == 'miR-6807-3p': continue
        if sMirName == 'miR-6838-5p': continue
        if sMirName == 'miR-7153-5p': continue

        list_cSeg = get_filtered_dotbracket_list (sInDir, sGenome, sChrID, list_cGFF3, sMirName, int(nWindowSize))

        # Text Output
        output_text(sOutDir, sMirName, list_cSeg)

        # Pickle Dump
        pickle.dump(list_cSeg, open('%s/%s.list' % (sOutDir, sMirName), 'wb'))

        print('Filtering %s %s...DONE' % (sMirName, sChrID))
    #loop END: sMirName

    print('Filtering Dot-Bracket...DONE %s' % time.ctime())
#def END: filter_dot_bracket


def get_filtered_dotbracket_list (sInDir, sGenome, sChrID, list_cGFF3, sMirName, nWindowSize):

    list_cSeg      = pickle.load(open('%s/%s/%s.list' % (sInDir, sChrID, sMirName), 'rb'))
    list_cSeg_filtered = []

    nKeep              = 0
    nDiscard           = 0

    for cSeg in list_cSeg:

        nSeedLength   = len(cSeg.sSeedSeq)

        # Check Dot-Bracket Notation : Unpaired Seed Site
        sSeed_dotbr   = cSeg.sDotBracket[cSeg.nSeedStart:cSeg.nSeedEnd]

        if sSeed_dotbr == '.' * nSeedLength:
            cSeg = check_seed_additional_pairing(cSeg, sGenome, list_cGFF3, sChrID[:-1], sMirName, nWindowSize)
            list_cSeg_filtered.append(cSeg)

            nKeep += 1

        else:
            nDiscard += 1
        #if END: sSeed_dotbr
    #loop END: cSeg

    print('Total Segments %d' % len(list_cSeg))
    print('Unpaired Seed %d'  % nKeep)
    print('Discarded %d'      % nDiscard)

    #V-S Check: Empty List
    if nKeep != len(list_cSeg_filtered):
        sys.exit('ERROR: list_sUnpairedSeed Size= %d' % len(list_cSeg_filtered))

    return list_cSeg_filtered
#def END: get_filtered_dotbracket_list


def check_seed_additional_pairing (cSeg, sGenome, list_cGFF3, sChrID, sMirName, nWindowSize):
    sGenomeDir = 'hg19' if sGenome == 'hsa' else 'mm9'
    sVertDir   = 'vertebrate_100way' if sGenome == 'hsa' else 'vertebrate'

    #cMiRNA Class Object with Mature Seq Info
    cMature         = get_miRNA_mature_seq (sGenome, sGenomeDir, sVertDir, sChrID, sMirName, list_cGFF3)

    sSeedRegion     = cSeg.sDotBracket[cSeg.nSeedStart:cSeg.nSeedEnd]
    sPairedRegion   = cSeg.sDotBracket[cSeg.nSeedStart-nMAX_INTERDISTANCE:cSeg.nSeedEnd]

    #V-S Check: Paired In Open:
    if '(' in sSeedRegion and ')' in sSeedRegion:
        sys.exit('ERROR: Invalid SeedSite Dot-Bracket= %s' % sSeedRegion)

    cSeg.nLowerWin = nWindowSize
    cSeg.nUpperWin = nWindowSize

    if ')' in sPairedRegion and '(' not in sPairedRegion:  # Find Dot-Bracket with Paired Before Seed


        ##### Check 8mer #####
        if cSeg.sSegmentSeq[cSeg.nSeedEnd] == 'A' and cSeg.sDotBracket[cSeg.nSeedEnd] == '.':
            cSeg           = copy.deepcopy(cSeg)
            cSeg.sSeedType = '8mer'
            cSeg.nSeedEnd  = cSeg.nSeedEnd + 1
            cSeg.sSeedSeq  = cSeg.sSeedSeq + 'A'
            cSeg.nLowerWin = nWindowSize
            cSeg.nUpperWin = nWindowSize - 1

        ##### Check Additional Seed #####
        sMaturSeq_afterseed = reverse_complement(cMature.sMatureSeq[8:]) # Get rev com mature sequence after seed
        bContinuous         = True

        for i in range(nADD_MATCH_MER_SIZE):

            nSeedStartIndex = cSeg.nSeedStart
            # Evaluate position 3' to 5' from site match
            sSegmentNuc = cSeg.sSegmentSeq[nSeedStartIndex-1-i] # Segment residue
            sDotNuc     = cSeg.sDotBracket[nSeedStartIndex-1-i] # DotBracket residue
            sMirNuc     = sMaturSeq_afterseed[-(i+1)]           # Mature miRNA residue

            if sSegmentNuc == sMirNuc and sDotNuc == '.':

                if bContinuous == True:
                    cSeg.nLowerWin  -= 1
                    nSeedStartIndex -= 1
                    cSeg.sSeedSeq    = sSegmentNuc + cSeg.sSeedSeq
                    cSeg.sSeedType   = '%smer' % len(cSeg.sSeedSeq)

            else:
                cSeg.nSeedStart = nSeedStartIndex
                bContinuous = False

            #if ENd: sSegmenNuc...
        #loop END: i

        ##### Recalculate Conservation Score #####

        # Genomic start and end positions
        nSeedStart = cSeg.nStartPos + cSeg.nLowerWin
        nSeedEnd   = cSeg.nStartPos + cSeg.nLowerWin + len(cSeg.sSeedSeq)

        list_fCons = fetch_cbc_cons(sGenomeDir, sVertDir, sChrID, nSeedStart, nSeedEnd)

        #V-S Check: Seed and list_fCons length
        if len(cSeg.sSeedSeq) != len(list_fCons):
            sys.exit('ERROR: Invalud list_fCons Size= %d\t Seed Length= %d' % (len(list_fCons), len(cSeg.sSeedSeq)))

        cSeg.fSeedCons  = sum(list_fCons)

        if bDEBUGMODE == 1:
            print(cSeg.sSeedType)
            print(' ' * cSeg.nLowerWin + cSeg.sSegmentSeq[cSeg.nSeedStart:cSeg.nSeedEnd])
            print(' ' * cSeg.nLowerWin + cSeg.sSeedSeq)
            print(' ' * (cSeg.nLowerWin - len(cMature.sMatureSeq) + len(cSeg.sSeedSeq)) + reverse_complement(cMature.sMatureSeq))
            print(cSeg.sSegmentSeq)
            print(cSeg.sDotBracket)
            print(' '* (cSeg.nSeedStart-nMAX_INTERDISTANCE) +  sPairedRegion)
    #if END: sPairedRegion...

    return cSeg
#def END:check_seed_additional_pairing


## endregion


## region Screen Dot-Bracket - Part 2: Apply Additional Constraints

def apply_additional_constraints (sInDir, sOutDir, sGenome, sChrID, nWindowSize):
    print('Applying Additional Constraints %s' % time.ctime())

    dict_sSeedSeq = load_87_conserved_mirna_seedseq()
    list_sMirName = sorted(list(dict_sSeedSeq.keys()))
    list_cGFF3    = parse_GFF3_files(sGenome)


    for sMirName in list_sMirName:

        if sMirName == 'miR-6766-3p': continue
        if sMirName == 'miR-6807-3p': continue
        if sMirName == 'miR-6838-5p': continue
        if sMirName == 'miR-7153-5p': continue

        print('Filtering %s' % sMirName)

        # Load Dot-Bracket List
        list_cSeg = pickle.load(open('%s/%s/%s.list' % (sInDir, sChrID, sMirName), 'rb'))

        # Constraint 1: Additional Downstream Unpaired Region
        list_cSeg = constraint1_additional_unpaired_region (list_cSeg, int(nWindowSize))


        # Constraint 2: Open Region Sequence Match with miRNA 3' Region
        list_cSeg = constraint2_match_with_miRNA_3prime_region (sGenome, sChrID[:-1], sMirName, list_cSeg, list_cGFF3, int(nWindowSize))

        print('list_cSeg', len(list_cSeg))

        # Pickle Dump
        pickle.dump(list_cSeg, open('%s/%s.list' % (sOutDir, sMirName), 'wb'))

        print('Filtering %s....DONE' % sMirName)
    #loop END: sMirName

    print('Applying Additional Constraints %s' % time.ctime())
#def END: filter_dot_bracket


def constraint1_additional_unpaired_region (list_cSeg, nWindowSize):

    list_sFiltered = []

    for cSeg in list_cSeg:

        nSeedLength     = len(cSeg.sSeedSeq)

        nSeedStartPos   = nWindowSize
        nSeedEndPos     = nWindowSize + nSeedLength

        sSeedRegion     = cSeg.sDotBracket[nSeedStartPos:nSeedEndPos]
        sPairedRegion   = cSeg.sDotBracket[nSeedStartPos-nMAX_INTERDISTANCE:nSeedStartPos]

        #V-S Check: Paired In Open:
        if '(' in sSeedRegion and ')' in sSeedRegion:
            sys.exit('ERROR: Invalid SeedSite Dot-Bracket= %s' % sSeedRegion)


        if ')' in sPairedRegion and '(' not in sPairedRegion:  # Find Dot-Bracket with Paired Before Seed

            sBeforeSeedRegion = cSeg.sDotBracket[:nSeedStartPos]    # Everything Before the Seed Site

            # Regular Expression : Find (....) Pattern with any number of .'s up to limit
            sSeedPattern      = '\(\.{%d,%d}?\.\)' % (nMIN_ADDITIONAL_SITE_SIZE, nMAX_ADDITIONAL_SITE_SIZE)

            list_sOpenSeq     = []

            for sReIndex in re.finditer(sSeedPattern, sBeforeSeedRegion):
                nIndexStart = int(sReIndex.start()) + 1  # '(' Not Included
                nIndexEnd   = int(sReIndex.end())   - 1  # ')' Not Included

                sOpenSeq    = cSeg.sSegmentSeq[nIndexStart:nIndexEnd]

                sKey        = '%s,%s,%s' % (nIndexStart, nIndexEnd, sOpenSeq)

                #if nIndexStart:

                if nSeedStartPos - nIndexEnd > nMAX_INTERDISTANCE: continue
                #print(sMiRNA, nSeedStartPos - nIndexEnd)

                list_sOpenSeq.append(sKey)
            #loop END: sReIndex

            if list_sOpenSeq:

                cSeg.list_sOpenSeq = list_sOpenSeq
                list_sFiltered.append(cSeg)
            #if END: list_sOpenSeq

        #if END: ')' in sPairedRegion and '(' not in sPairedRegion
    #loop END: cSeg

    return list_sFiltered
#def END: additional_unpaired_region


def constraint2_match_with_miRNA_3prime_region (sGenome, sChrID, sMirName, list_cSeg, list_cGFF3, nWindowSize):

    sGenomeDir = 'hg19' if sGenome == 'hsa' else 'mm9'
    sVertDir   = 'vertebrate_100way' if sGenome == 'hsa' else 'vertebrate'

    #cMiRNA Class Object with Mature Seq Info
    cMature        = get_miRNA_mature_seq (sGenome, sGenomeDir, sVertDir, sChrID, sMirName, list_cGFF3)

    dict_sMature   = {}

    dict_sMature   = get_motif_dictionary (dict_sMature, cMature.sMatureSeq[8:])

    list_sFiltered = []

    for cSeg in list_cSeg:
        cSeg = copy.deepcopy(cSeg)
        nSeedLength      = len(cSeg.sSeedSeq)

        dict_sOpenRegion = {}

        for sKey in cSeg.list_sOpenSeq:

            nIndexStart      = int(sKey.split(',')[0])

            sOpenSeq         = reverse_complement(sKey.split(',')[2])

            dict_sOpenRegion = get_motif_dictionary(dict_sOpenRegion, sOpenSeq)

            for sOpenRegion in dict_sOpenRegion:
                dict_sOpenRegion[sOpenRegion][0] = dict_sOpenRegion[sOpenRegion][0] + nIndexStart
                dict_sOpenRegion[sOpenRegion][1] = dict_sOpenRegion[sOpenRegion][1] + nIndexStart
            #loop END: sOpenRegion
        #loop END: sKey

        # Find Matching Sequences
        list_sAddPairSeq = list(dict_sMature.keys() & dict_sOpenRegion.keys())
        list_sAddPairSeq = sorted(list_sAddPairSeq, key=len, reverse=True)


        if not list_sAddPairSeq: continue # if additional pairing

        list_sAdd_Cons = []

        for sAddPairSeq in list_sAddPairSeq[:1]:

            nIndexStart, nIndexEnd = dict_sOpenRegion[sAddPairSeq]

            if cSeg.sStrand == '+':
                nStart = cSeg.nStartPos + nIndexStart
                nEnd   = cSeg.nStartPos + nIndexEnd
            else:

                nStart = cSeg.nEndPos - nIndexEnd
                nEnd   = cSeg.nEndPos - nIndexStart
            #if END: sChrID

            list_fCons          = fetch_cbc_cons (sGenomeDir, sVertDir, cSeg.sChrID, nStart-1, nEnd)


            #V-S Check: List Size
            if len(list_fCons) != len(sAddPairSeq):
               sys.exit('ERROR: Invalid Conservation List Size= %d\tStart= %d\tEnd= %d\tAddPairSeq= %s\t Cons= %s' % (len(list_fCons), nStart, nEnd, sAddPairSeq, sum(list_fCons)))

            cSeg.fConsScore     = cSeg.fSeedCons + sum(list_fCons)

            cSeg.sAddPair       = sAddPairSeq

            if cSeg.fConsScore < nSECOND_CONS_FILTER: continue  ## Conservation Score Cut-off

            list_sAdd_Cons.append([sAddPairSeq, cSeg.fConsScore])

        #loop END: sMatchedSeed

        if list_sAdd_Cons:
            cSeg.list_sAdd_Cons = list_sAdd_Cons
            list_sFiltered.append(cSeg)
        #if END: list_sMatchedSeed

    #loop END: sChrID , sMiRNA, sSeedSite, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket, list_sOpenSeq

    return list_sFiltered
#def pairing_with_miRNA_3prime_region


def constraint2_match_with_miRNA_3prime_region_v2 (sGenome, cChrID, sMirName, list_cSeg, list_cGFF3, nWindowSize):

    sGenomeDir = 'hg19' if sGenome == 'hsa' else 'mm9'
    sVertDir   = 'vertebrate_100way' if sGenome == 'hsa' else 'vertebrate'

    #cMiRNA Class Object with Mature Seq Info
    cMature        = get_miRNA_mature_seq (sGenome, sGenomeDir, sVertDir, cChrID, sMirName, list_cGFF3)

    dict_sMature   = {}

    sMaskedSeq     = ('*' * 8) + cMature.sMatureSeq[8:]

    dict_sMature   = get_motif_dictionary(dict_sMature, sMaskedSeq, nADD_MATCH_MER_SIZE)

    list_sFiltered = []

    for cSeg in list_cSeg:

        nSeedLength      = len(cSeg.sSeedSeq)

        dict_sOpenRegion = {}

        for sKey in cSeg.list_sOpenSeq:

            nIndexStart      = int(sKey.split(',')[0])
            sOpenSeq         = reverse_complement(sKey.split(',')[2])

            dict_sOpenRegion = get_motif_dictionary(dict_sOpenRegion, sOpenSeq, nADD_MATCH_MER_SIZE)

            for sOpenRegion in dict_sOpenRegion:
                dict_sOpenRegion[sOpenRegion][0] = dict_sOpenRegion[sOpenRegion][0] + nIndexStart
                dict_sOpenRegion[sOpenRegion][1] = dict_sOpenRegion[sOpenRegion][1] + nIndexStart
            #loop END: sOpenRegion
        #loop END: sKey

        # Find Matching Sequences
        list_sMatchedSeed = list(dict_sMature.keys() & dict_sOpenRegion.keys())

        if not list_sMatchedSeed: continue # if additional pairing

        list_sAdd_Cons = []

        for sMatchedSeed in list_sMatchedSeed:

            nIndexStart, nIndexEnd = dict_sOpenRegion[sMatchedSeed]


            if cSeg.sStrand == '+':
                nStart = cSeg.nStartPos + nIndexStart
                nEnd   = cSeg.nStartPos + nIndexStart
            else:
                nStart = cSeg.nEndPos - nIndexEnd
                nEnd   = cSeg.nEndPos - nIndexStart
            #if END: sChrID

            cSeg.list_sCons_add = fetch_cbc_cons (sGenomeDir, sVertDir, cSeg.sChrID, nStart-1, nEnd)

            cSeg.fConsScore     = cSeg.fSeedCons + sum(cSeg.list_sCons_add)

            if cSeg.fConsScore < nSECOND_CONS_FILTER: continue  ## Conservation Score Cut-off

            list_sAdd_Cons.append([sMatchedSeed, cSeg.fConsScore])
        #loop END: sMatchedSeed

        if list_sAdd_Cons:
            cSeg.list_sAdd_Cons = list_sAdd_Cons
            list_sFiltered.append(cSeg)
        #if END: list_sMatchedSeed

    #loop END: sChrID , sMiRNA, sSeedSite, nPos, fMFE, fLocalAU, sSegmentSeq, sDotBracket, list_sOpenSeq

    return list_sFiltered
#def pairing_with_miRNA_3prime_region


def get_miRNA_mature_seq (sGenome, sGenomeDir, sVertDir, sChrID, sMirName, list_cGFF3):

    try:
        cMatureMir            = { cMir.sMirName.lower().replace('%s-' % sGenome, ''): cMir for cMir in list_cGFF3
                                  if cMir.sFeature == 'miRNA'}[sMirName.lower()]
    except KeyError:

        if   sMirName == 'miR-203a':                 sMirName = 'mir-203a-5p'
        elif sMirName == 'miR-219a-5p':              sMirName = 'mir-219-5p'
        elif sMirName == 'miR-429':                  sMirName = 'miR-429-5p'
        elif '-5p' in sMirName or '-3p' in sMirName: sMirName = sMirName[:-3]

        cMatureMir            = { cMir.sMirName.lower().replace('%s-' % sGenome, ''): cMir for cMir in list_cGFF3
                                   if cMir.sFeature == 'miRNA'}[sMirName.lower()]

    cGenome               = cFasta(dict_sGENOME_DIR[sGenome])


    cMatureMir.sMatureSeq = cGenome.fetch(cMatureMir.sChrID, cMatureMir.nStartPos-1,
                                          cMatureMir.nEndPos, cMatureMir.sStrand)


    cMatureMir.list_fCBC  = fetch_cbc_cons (sGenomeDir, sVertDir, cMatureMir.sChrID,
                                            cMatureMir.nStartPos-1, cMatureMir.nEndPos)

    return cMatureMir
#def END: get_miRNA_mature_seq


def get_motif_dictionary (dict_sSeq, sSeq):

    for i in range(len(sSeq)):

        for j in range(len(sSeq) - nADD_MATCH_MER_SIZE + 1):

            sKeySeq = sSeq[i:j + nADD_MATCH_MER_SIZE]

            if len(sKeySeq) < nADD_MATCH_MER_SIZE: continue

            if sKeySeq not in dict_sSeq:
                dict_sSeq[sKeySeq] = []

            dict_sSeq[sKeySeq] = [len(sSeq) - (j+nADD_MATCH_MER_SIZE), len(sSeq) - 1 - i]
            #dict_sSeq[sKeySeq] = [i, (j+nADD_MATCH_MER_SIZE)]
        #loop END: j
    #loop END: i

    #V-S Check: Empty Dictionary
    if not dict_sSeq:
        sys.exit('ERROR: dict_sSeq Size= %d' % len(dict_sSeq))

    return dict_sSeq
#def END: get_motif_dictionary


def get_motif_dictionary_v2 (dict_sSeq, sSeq, sSeqLen):

    for i in range(len(sSeq) - sSeqLen + 1):

        sKeySeq = sSeq[i:(i+sSeqLen)]

        if '*' in sKeySeq: continue # SKIP Masked Region

        if sKeySeq not in dict_sSeq:
            dict_sSeq[sKeySeq] = []

        dict_sSeq[sKeySeq] = [i, i+sSeqLen]
    #loop END: i

    #V-S Check: Empty Dictionary
    if not dict_sSeq:
        sys.exit('ERROR: dict_sSeq Size= %d' % len(dict_sSeq))

    return dict_sSeq
#def END: get_motif_dictionary

## endregion


## region Sort and Evaluate - RNAFold and RNAPlot

def sort_and_evaluate_segments (sInDir, sOutDir, sGenome,  sMirName, nWindowSize):
    print('Sorting and Evaluating Segments %s' % time.ctime())

    list_cSeg = get_compiled_dotbracket_list(sInDir, sGenome, sMirName)
    #list_cSeg    = classify_list (list_cSeg)

    # Pickle Dump
    pickle.dump(list_cSeg, open('%s/%s.class' % (sOutDir, sMirName), 'wb'))

    print('Sorting and Evaluating Segments...DONE %s' % time.ctime())
#def sort_and_evaluate_segments


def get_compiled_dotbracket_list (sInDir, sGenome,  sMirName):

    cGenome = cFasta(dict_sGENOME_DIR[sGenome])

    list_sChrIDs       = ['%s%s' % (sChrID, sStrand) for sChrID in cGenome.sChrIDList for sStrand in ['+', '-']]
    list_cSeg      = []

    for sChrID in list_sChrIDs:
        if 'M' in sChrID: continue   # SKIP ChrM
        list_cSeg += pickle.load(open('%s/%s/%s.list' % (sInDir, sChrID, sMirName), 'rb'))
    #loop END: sChrID

    return sorted(list_cSeg, key=lambda c:c.fSeedCons, reverse=True)
#def END: get_compiled_dotbracket_list


def classify_list (list_sDotBracket):
    print('Classifying List')
    list_cSeg = []

    for sSegmentData in list_sDotBracket:

        cSeg = cSegment()

        cSegment_assign_segment_data(cSeg, sSegmentData)

        list_cSeg.append(cSeg)
    #loop END: sSegmentData

    print('Classifying List...DONE')
    return sorted(list_cSeg, key=lambda c:(c.fMFE, -c.fLocalAU))
#def END: classify_list

## endregion


## region Compile Sorted List and Plot SS

def full_compile_list (sInDir, sOutDir, nWindowSize, sTargetSeedType, sSortBy):
    print('Compiling Full List %s' % time.ctime())

    dict_sSeedSeq = load_87_conserved_mirna_seedseq()
    list_sMirName = sorted(list(dict_sSeedSeq.keys()))

    list_cSeg     = []

    for sMirName in list_sMirName:
        if sMirName == 'miR-6766-3p': continue
        if sMirName == 'miR-6807-3p': continue
        if sMirName == 'miR-6838-5p': continue
        if sMirName == 'miR-7153-5p': continue


        list_cSeg += pickle.load(open('%s/%s.class' % (sInDir, sMirName), 'rb'))

    #loop END: sMirName

    if sSortBy == 'Cons':
        list_cSeg = sorted(list_cSeg, key=lambda c:c.fConsScore, reverse=True)

    elif sSortBy == 'PairLength':
        list_cSeg = sorted(list_cSeg, key=lambda c:((len(c.sSeedSeq) + len(c.sAddPair), c.fConsScore)), reverse=True)


    # Plot Top N Secondary Structure on SVG
    list_cSeg = plot_SS_on_SVG_individual (sOutDir, list_cSeg, 'SVG_%scons_%sopen_%s'
                % (nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE, sSortBy), int(nWindowSize), sTargetSeedType, sSortBy)

    # Pickle Dump
    pickle.dump(list_cSeg, open('%s/%s_%scons_%sopen_%s.class'
                                    % (sOutDir, sTargetSeedType, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE, sSortBy), 'wb'))

    print('Compiling Full List...DONE %s' % time.ctime())
#def END: full_compile_list

## endregion


## region Extend Segment and Plot SS

def extended_ss (sInDir, sOutDir, nWindowSize, sTargetSeedType, sSortBy):
    print('Compiling Full List %s' % time.ctime())

    list_cSeg = pickle.load(open('%s/%s_%scons_%sopen_%s.class' % (sInDir, sTargetSeedType, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE, sSortBy), 'rb'))

    if sSortBy == 'Cons':
        list_cSeg = sorted(list_cSeg, key=lambda c:c.fConsScore, reverse=True)

    elif sSortBy == 'PairLength':
        list_cSeg = sorted(list_cSeg, key=lambda c:(len(c.sSeedSeq) + len(c.sAddPair)), reverse=True)


    # Plot Top N Secondary Structure on SVG
    plot_SS_on_SVG_extended (sOutDir, list_cSeg, 'SVG_%scons_%sopen_%s_ext'
                % (nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE, sSortBy), int(nWindowSize), sTargetSeedType, sSortBy)

    print('Compiling Full List...DONE %s' % time.ctime())
#def END: full_compile_list


def extend_segment_info (cSeg, nWindowSize, nEXTEND_LEN, sTargetSeedType):

    cSeg            = copy.deepcopy(cSeg)

    cGenome         = cFasta(dict_sGENOME_DIR['hsa'])

    nSeedLength     = len(cSeg.sSeedSeq)


    nStartPos       = int(cSeg.nStartPos) - (nEXTEND_LEN - 1)
    nEndPos         = int(cSeg.nEndPos)   + (nEXTEND_LEN)

    cSeg.sSegmentSeq_ext = cGenome.fetch(cSeg.sChrID, nStartPos-1, nEndPos, cSeg.sStrand)

    cSeg.sDotBracket_ext, cSeg.fMFE_ext = get_dotbracket_and_MFE (cSeg.sSegmentSeq_ext, 0.0)

    nSeedStartPos        = cSeg.nLowerWin + nEXTEND_LEN
    nSeedEndPos          = cSeg.nUpperWin + nEXTEND_LEN + nSeedLength


    if '(' in cSeg.sDotBracket_ext[nSeedStartPos:nSeedEndPos] or ')' in cSeg.sDotBracket_ext[nSeedStartPos:nSeedEndPos]:
        bOpenSeedFlag = False

    else: bOpenSeedFlag = True

    # Get LocalAU Score : Code Credit - Doyeon Kim
    cSeg.fLocalAU_ext        = get_LocalAU_score (nSeedStartPos, nSeedEndPos, cSeg.sSegmentSeq_ext,
                                                  nSeedLength, nWindowSize)

    return cSeg, bOpenSeedFlag
#def END: extend_segment_info

## endregion


## region Overlapped
def compare_lists (sInDir, sOutDir):

    print('%s/75bps/full_compile_list/full_compile_list_%scons_%sopen.class' %(sInDir, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE))

    list_sSegment_75bps = pickle.load(open('%s/75bps/full_compile_list/full_compile_list_%scons_%sopen.class'  % (sInDir, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE), 'rb'))

    list_sSegment_100bps = pickle.load(open('%s/100bps/full_compile_list/full_compile_list_%scons_%sopen.class'% (sInDir, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE), 'rb'))

    list_sSegment_125bps = pickle.load(open('%s/125bps/full_compile_list/full_compile_list_%scons_%sopen.class'% (sInDir, nSECOND_CONS_FILTER, nADD_MATCH_MER_SIZE), 'rb'))

    print('75bps', len(list_sSegment_75bps))
    print('100bps', len(list_sSegment_100bps))
    print('125bps', len(list_sSegment_125bps))

    list_sShared  = get_shared_list(list_sSegment_75bps, list_sSegment_100bps, list_sSegment_125bps)

    list_cSeg = [cSeg for cSeg in list_sSegment_75bps if '%s:%s-%s' % (cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd) in list_sShared]

    print(len(list_sShared))
    print(len(list_cSeg))

    for cSeg in list_cSeg:
        sMatchedSeq, fConserved = cSeg.list_sCons[0]

        if fConserved > 10:
            print(cSeg.sChrID+cSeg.sStrand, cSeg.sMirName, cSeg.fMFE, cSeg.fLocalAU, fConserved,  cSeg.nStartPos, cSeg.nEndPos, '%s:%s-%s' % (cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd))

            print(cSeg.sSegmentSeq)
    sys.exit()

    pickle.dump(list_cSeg, open('%s/%s_%scons_%sopen.class' % (sOutDir, 'full_overlapped_list',
                                                                  nMIN_CONSERVE_SCORE, nADD_MATCH_MER_SIZE), 'wb'))

    plot_SS_on_SVG_v2 (sOutDir, list_cSeg, 'full_compiled_list_%scons_%sopen_overlapped'
                       % (nMIN_CONSERVE_SCORE, nADD_MATCH_MER_SIZE), 75, 'overlap')
#def END: get_common_segments


def get_shared_list (list_sSegment_75bps, list_sSegment_100bps, list_sSegment_125bps):

    list_75bps = ['%s:%s-%s' % (cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd) for cSeg in list_sSegment_75bps]

    list_100bps = ['%s:%s-%s' % (cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd) for cSeg in list_sSegment_100bps]

    list_125bps = ['%s:%s-%s' % (cSeg.sChrID, cSeg.nWindowStart, cSeg.nWindowEnd) for cSeg in list_sSegment_125bps]

    return list(set(list_75bps) & set(list_100bps) & set(list_125bps))
#def END: get_shared_list

## endregion


## region Expression Profile

# sRNA Seq - Human Cell-Line
def get_maturemir_RPM (sAnalysis, sCellLineID, sInDir, sOutDir):
    print('Getting sRNA RPM List %s' % (time.ctime()))

    sTartgetMir     = 'mir-302a-3p'

    list_cGFF3      = parse_GFF3_files(sAnalysis)

    sInputFile_aln  = '%s/%s.aln'          % (sInDir, sCellLineID)


    list_sOut       = []

    for cMir in list_cGFF3:

        if cMir.sFeature != 'miRNA_primary_transcript': continue
        list_cGFF3_mature = [cGFF3 for cGFF3 in list_cGFF3 if cGFF3.sFeature != 'miRNA_primary_transcript'
                             and include(cMir.sChrID, cGFF3.sChrID, [cGFF3.nStartPos, cGFF3.nEndPos], [cMir.nStartPos, cMir.nEndPos]) ==  0 ]

        for cMir_mature in list_cGFF3_mature:

            nReadDataList   = access_aln_file(sInputFile_aln, cMir_mature.sChrID, cMir_mature.sStrand, cMir_mature.nStartPos, cMir_mature.nEndPos, 1, 'expr')

            fRPM            = sum([ nUniqCnt + nMultiCnt for nStartPos, nEndPos, nReadCnt, nUniqCnt, nMultiCnt in nReadDataList ])

            list_sOut.append([cMir_mature.sMirName, fRPM])


            if sTartgetMir in cMir_mature.sMirName:

                pickle.dump([cMir_mature.sMirName, fRPM], open('%s/%s_%s_norm.rpm' % (sOutDir, sCellLineID, sTartgetMir), 'wb'))

        #loop END: cMir_mature
    #loop END: cMiRNA

    list_sOut = sorted(list_sOut, key=lambda e:e[1], reverse=True)

    OutFile         = open('%s/%s_RPM.txt' % (sOutDir, sCellLineID), 'w')

    for sMatureMir, fRPM in list_sOut:
        sOut            = '%s\t%s\n' % (sMatureMir, fRPM)

        OutFile.write(sOut)
    #loop END: sMatureMir, fRPM
    OutFile.close()

    print('Getting sRNA RPM List...DONE %s' % (time.ctime()))
#def END: get_maturemir_RPM


def compile_targetmir_rpm_list (sAnalysis, sInDir, sOutDir):

    sTargetMir      = 'mir-302a-3p'

    list_sFiles      = get_file_id_list(sAnalysis, 'sRNA')

    dict_fRPM_hsa    = get_RPM_dictionary(sAnalysis, sInDir, sTargetMir, 'hsa', list_sFiles)
    dict_fRPM_mmu    = get_RPM_dictionary(sAnalysis, sInDir, sTargetMir, 'mmu', list_sFiles)

    list_sOutput_hsa = get_RPM_outlist (dict_fRPM_hsa)

    list_sOutput_mmu = get_RPM_outlist (dict_fRPM_mmu)


    print('*******************Human*******************')
    for sCellLine, fAvgRPM, fAvgPerRank, fAvgTotRPM in list_sOutput_hsa:
        print('%s\t%s\t%0.1f%%\t%s' % (sCellLine, fAvgRPM, fAvgPerRank, fAvgTotRPM))
    #loop END: sFileID, sCellLine, fRPM
    get_basic_stats (list_sOutput_hsa)

    print('*******************Mouse*******************')
    for sCellLine, fAvgRPM, fAvgPerRank, fAvgTotRPM in list_sOutput_mmu:
        print('%s\t%s\t%0.1f%%\t%s' % (sCellLine, fAvgRPM, fAvgPerRank, fAvgTotRPM))
    #loop END: sFileID, sCellLine, fRPM
    get_basic_stats (list_sOutput_mmu)
#def END: compile_targetmir_rpm_list


def get_RPM_dictionary (sAnalysis, sInDir, sTargetMir, sTarGenome, list_sFiles):

    dict_fRPM        = {}

    for sFileName in list_sFiles:
        sMirName, fRPM  = pickle.load(open('%s/%s_%s_norm.rpm' % (sInDir, sFileName, sTargetMir), 'rb'))


        fAvgRPM_allmirs, fPercentRank = get_all_mir_stats(sInDir, sFileName, sTargetMir)

        sGenome         = assign_cell_line_info(sFileName)['genome']
        sCellLine       = assign_cell_line_info(sFileName)['cell_line']

        if sGenome == sTarGenome:
            if sCellLine not in dict_fRPM:
                dict_fRPM[sCellLine] = []

            dict_fRPM[sCellLine].append([fRPM, fPercentRank, fAvgRPM_allmirs])
    #loop END: sFileName

    #V-S Check: Empty Dictionary
    if not dict_fRPM:
        sys.exit('ERROR: dict_fRPM Size= %d' % len(dict_fRPM))

    return dict_fRPM
#def END: get_RPM_dictionary


def get_all_mir_stats (sInDir, sFileName, sTargetMir):

    sInFileDir   = '%s/%s_RPM.txt' % (sInDir, sFileName)

    nCnt         = sum([1 for i in open(sInFileDir)])

    dict_fRPM    = {}

    fPercentRank = 0.0

    InFile       = open(sInFileDir)

    for i, sReadLine in enumerate(InFile):
        #File Format:
        #Column:        1               |  2
        #ID:            sMirName        |  fRPM
        #Example:       hsa-mir-6723-5p |  3.66542704025904

        list_sColumn = sReadLine.strip('\n').split('\t')

        sMirName = list_sColumn[0]
        fRPM     = float(list_sColumn[1])

        if sTargetMir in sMirName:
            fPercentRank = (1 - (i/nCnt)) * 100

        if sMirName not in dict_fRPM:
            dict_fRPM[sMirName] = 0.0

        dict_fRPM[sMirName] = fRPM
    #loop END: sMirName, fRPM
    InFile.close()

    #V-S Check: No PercentRank
    if fPercentRank == 0:
        sys.exit('ERROR: fPercentRank= %s\t sFileName= %s' % (fPercentRank, sFileName))

    return sum([dict_fRPM[sMirName] for sMirName in dict_fRPM]) / len(dict_fRPM), fPercentRank
#def END: get_average_rpm


def get_basic_stats (list_sOutput):

    #V-S Check: List Size = 0
    if not list_sOutput:
        sys.exit('ERROR: list_sOutput Size= %d' % len(list_sOutput))

    fAvgRPM = sum([fRPM for sCellLine, fRPM, fPerRank, fAvgTotRPM in list_sOutput]) / len(list_sOutput)
    fMaxRPM = max([fRPM for sCellLine, fRPM, fPerRank, fAvgTotRPM in list_sOutput])
    fMinRPM = min([fRPM for sCellLine, fRPM, fPerRank, fAvgTotRPM in list_sOutput])

    print('Average: %s' % fAvgRPM)
    print('Max: %s' % fMaxRPM)
    print('Min: %s' % fMinRPM)
#def END: get_basic_stats


def get_RPM_outlist (dict_fRPM):

    list_sOutput = []

    for sCellLine in dict_fRPM:

        fAvgRPM     = sum([fRPM for fRPM, fPercentRank, fAvgRPM_allmirs in dict_fRPM[sCellLine]]) / len(dict_fRPM[sCellLine])
        fAvgPerRank = sum([fPercentRank for fRPM, fPercentRank, fAvgRPM_allmirs in dict_fRPM[sCellLine]]) / len(dict_fRPM[sCellLine])
        fAvgTotRPM  = sum([fAvgRPM_allmirs for fRPM, fPercentRank, fAvgRPM_allmirs in dict_fRPM[sCellLine]]) / len(dict_fRPM[sCellLine])
        list_sOutput.append([sCellLine, fAvgRPM, fAvgPerRank, fAvgTotRPM])
    #loop END: sCellLine

    list_sOutput = sorted(list_sOutput, key=lambda e:e[1], reverse=True)

    #V-S Check: Empty List

    if not list_sOutput:
        sys.exit('ERROR: list_sOutput Size= %d' % len(list_sOutput))

    return list_sOutput
#def END: get_RPM_outlist


def get_target_mir_rpkm_list (sGenome, sAnalysis, sStudy, sInDir, sOutDir):

    sFileIDName = '%s_%s_%s' % (sGenome, sAnalysis, sStudy)

    list_sFiles     = get_file_id_list(sFileIDName, 'mRNA')

    dict_fRPKM      = {}

    list_sTargetGeneSym   = ['SEGMENT302o', 'SEGMENT302_100', 'SEGMENT302_250', 'SEGMENT302_500', 'SEGMENT302_1000',
                             'SEGMENT_NEGCON', 'EEEF1A1', 'GAPDH', 'ACTB']

    for sFileName in list_sFiles:

        if sFileName not in dict_fRPKM:
            dict_fRPKM[sFileName] = []

        sInFileDir = '%s/%s/accepted_hits.sorted.rpkm' % (sInDir, sFileName)

        dict_fRPKM[sFileName] = load_rpkm_list(sInFileDir, list_sTargetGeneSym)
    #loop END: sFileName

    list_sOut           = []

    for sFileName in dict_fRPKM:

        list_sRPKM = '\t'.join([('%0.3f' % dict_fRPKM[sFileName][sGeneSym][0]) for sGeneSym in list_sTargetGeneSym])

        list_sOut.append([sFileName, list_sRPKM])
    #loop END: sFileName

    list_sOut = sorted(list_sOut, key=lambda e:e[0])

    sHeader = 'GENOME\tFILE_ID\tCELLSTAGE\t' + '\t'.join(list_sTargetGeneSym)
    print(sHeader)

    for sFileName, list_sRPKM in list_sOut:

        sGenome   = assign_cell_line_info(sFileName)['genome']
        sCellLine = assign_cell_line_info(sFileName)['cell_line']

        sOut = '%s\t%s\t%s\t%s' % (sGenome, sFileName, sCellLine, list_sRPKM)
        print(sOut)
    #loop END: sFileName, list_sRPKM
#def END:  get_target_mir_rpkm_list


def load_rpkm_list (sInFileDir, list_sTargetGene):

    nCnt          = sum([1 for i in open(sInFileDir)])

    dict_sSegment = {}

    InFile        = open(sInFileDir)
    for i, sReadLine in enumerate(InFile):
        # RPKM File Format:
        # Column Number:        0          | 1            | 2
        # Column Description:   GeneSym    | NMID         | RPKM
        # Column Example:       SFTPA2     | NM_001098668 | 14431.869354

        list_sColumn = sReadLine.strip('\n').split('\t')

        #V-S Check: Number of Columns
        if len(list_sColumn) != 3:
            sys.exit('ERROR: list_sColumn Size= %d\t Columns= %s' % (len(list_sColumn), list_sColumn))

        sGeneSym = list_sColumn[0].upper()
        sNMID    = list_sColumn[1]
        fRPKM    = float(list_sColumn[2])

        if sGeneSym in list_sTargetGene:

            if sGeneSym not in dict_sSegment:
                dict_sSegment[sGeneSym] = []

            dict_sSegment[sGeneSym] = [fRPKM, (1 - i/nCnt) * 100]
    #loop END: sReadLine
    InFile.close()

    #V-S Check: Dictionary Size
    if len(dict_sSegment) != len(list_sTargetGene):
        sys.exit('ERROR: Invalid Found in %s\n dict_sSegment Size= %d list_sTartgetGeneSym Size= %d' %
                 (sInFileDir, len(dict_sSegment), len(list_sTargetGene)))

    return dict_sSegment
#def END: load_rpkm_list

## endregion



def main():
    print('This is the main')
#def END: main()

if __name__ == '__main__':
    if len(sys.argv) == 1: main()
    else:
        function_name       = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys(): locals()[function_name](*function_parameters)
        else: sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    #if END: len(sys.argv)
#if END: __name__
