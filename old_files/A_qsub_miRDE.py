#!/extdata6/Sukjun/opt/python3/bin/python3

bDEBUGMODE = 0
bEXTENDED  = 1
import os, sys, itertools, pickle, time, array, re, copy
from collections import Counter
from scipy import stats
from numpy import array, empty

## region Globals

sBASE_DIR         = '/extdata5/Jinman/05_DecayElement'
sALN_DIR          = '/extdata5/Jinman/00_CG_srSeq/04_Core_Analysis/04_indexed_sam'
sFILEID_DIR_sRNA  = '/extdata5/Jinman/00_CG_srSeq/05_File_IDs'
sFILEID_DIR_mRNA  = '/extdata5/Jinman/03_CG_mrSeq'
fPSEUDO_COUNT     = float('1.0e-300')


list_sCHRID_HSA   = [   '1','2','3','4','5','6','7','8','9','10',
                        '11','12','13','14','15','16','17','18','19','20',
                        '21', '22','X','Y']

dict_sGENOME_DIR  = {   'hsa':'/extdata5/Jinman/cancermir/genome/hg19/hg19.fa',
                        'mmu':'/extdata5/Jinman/reference_genome/mm10/mm10.fa',
                        'mml':'/extdata5/Jinman/reference_genome/MMUL1.0/MMUL1.0.fa',
                        'rno':'/extdata5/Jinman/reference_genome/rn4/rn4.fa',
                        'bta':'/extdata5/Jinman/reference_genome/bosTau6/bosTau6.fa',
                        'gal':'/extdata5/Jinman/reference_genome/galGal4/galGal4.fa',
                        'TCGA':'/extdata5/Jinman/cancermir/genome/hg19/hg19.fa'}

dict_sMIRBASE_DIR = {   'hsa': '/extdata5/Jinman/reference_genome/mirbase/hsa.liftover.gff3',
                        'mmu': '/extdata5/Jinman/reference_genome/mirbase/mmu.liftover.gff3',
                        'mml': '/extdata5/Jinman/reference_genome/mirbase/mml.liftover.gff3',
                        'bta': '/extdata5/Jinman/reference_genome/mirbase/bta.liftover.gff3',
                        'rno': '/extdata5/Jinman/reference_genome/mirbase/rno.liftover.gff3',
                        'gal': '/extdata5/Jinman/reference_genome/mirbase/gga.liftover.gff3'}

nCONSERVED_MIRNA_COLUMN_SIZE = 2
nCONSERVED_MIRNAS            = 87
nSEED_FLANK_LENGTH           = 50

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
def copy_temp_core_script ():
    os.system('cp %s/B_core_miRDE.py %s/temp/B_core_miRDE_%s.py' % (sBASE_DIR, sBASE_DIR, sTIME_STAMP))
#def END: copy_temp_core_script


def reverse_complement(sInputSeq):

    dict_CompBases = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}

    list_sSeq      = list(sInputSeq) # Turns the sequence in to a gigantic list

    list_sComSeq   = [dict_CompBases[sBase] for sBase in list_sSeq]
                     # For each base, the dictionary value for that base key is assigned

    return ''.join(list_sComSeq)[::-1] # Empty string, join the list of bases, [start : end : backwards]
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


def get_file_id_list(sAnalysis, bSeqType):
    if bSeqType == 'sRNA': sInDir = sFILEID_DIR_sRNA
    else: sInDir = sFILEID_DIR_mRNA
    return [id[:-1] for id in open('%s/fileid_%s.txt' %(sInDir, sAnalysis), 'r')]
#def END: get_file_id_list


def get_dict_files (list_sFiles, bTest):
    dict_sFiles = {}

    for sFileID in list_sFiles:
        if bTest == False:
            sKey    = sFileID.split('-')[2]  # PatID  ex) 6743
        else:
            sKey    = sFileID.split('_')[0]  # mir-155 and Mock

        if sKey not in dict_sFiles:
            dict_sFiles[sKey] = None
        dict_sFiles[sKey] = sFileID     # Barcode
    #loop END: sFileID


    return dict_sFiles
#def END: get_dict_files


def load_87_conserved_mirna_seedseq ():
    print('Loading 87 Conversed miRNA Seed Sequence')
    dict_sSeedSeq = {}
    InFile = open('%s/87_conserved_mirna_7m8.txt' % sBASE_DIR, 'r')

    for sReadLine in InFile:
        # File Format
        # Column Number:     | 0          | 1       |
        # Column Descriptio: | miRNA      | 7mer-m8 |
        # Column Example:    | hsa-let-7a | GAGGTAG |

        list_sColumn = sReadLine.strip('\n').split('\t')

        # V-S Check:
        if len(list_sColumn) != nCONSERVED_MIRNA_COLUMN_SIZE:
            sys.exit('ERROR: Conserved MiRNA Column Size= %d\tContent= %s' % (len(list_sColumn), list_sColumn))

        sMirName    = list_sColumn[0]
        s7m8Seq     = reverse_complement(list_sColumn[1])
        s7A1Seq     = s7m8Seq[1:] + 'A'
        s8merSeq    = s7m8Seq + 'A'
        s6merSeq    = s7m8Seq[1:]

        if sMirName not in dict_sSeedSeq:
            dict_sSeedSeq[sMirName] = {}

        dict_sSeedSeq[sMirName] = {'8mer':s8merSeq, '7m8':s7m8Seq, '7A1':s7A1Seq, '6mer':s6merSeq}
    #loop END: sReadLine

    #V-S Check: Dictionary Size
    if len(dict_sSeedSeq) != nCONSERVED_MIRNAS:
        sys,exit('ERROR: dict_sSeedSeq Size= %d' % len(dict_sSeedSeq))

    print('Loading 87 Conversed miRNA Seed Sequence...DONE')

    return dict_sSeedSeq
#def END: load_87_conserved_mirna_seedseq


## endregion

## region class cFasta
re_nonchr = re.compile('[^a-zA-Z]')
class cFasta:
    def __init__(self, sRefFile):

        #V-S Check
        if not os.path.isfile(sRefFile):
           sys.exit('(): File does not exist')

        self.InFile     = open(sRefFile, 'r')
        self.sChrIDList = []
        self.nChromLen  = []
        self.nSeekPos   = []
        self.nLen1      = []
        self.nLen2      = []

        #V-S Check
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
        print('Fetching Chromosome Sequence %s' % sChrom)

        assert sChrom in self.sChrIDList, sChrom
        nChrom = self.sChrIDList.index(sChrom)

        if nFrom == None: nFrom = 0
        if nTo   == None: nTo = self.nChromLen[nChrom]

        assert(0 <= nFrom) and(nFrom < nTo) and(nTo <= self.nChromLen[nChrom])

        nBlank = self.nLen2[nChrom] - self.nLen1[nChrom]

        nFrom  = int(nFrom +(nFrom / self.nLen1[nChrom]) * nBlank) # Start Fetch Position

        nTo    = int(nTo   +(nTo   / self.nLen1[nChrom]) * nBlank) # End Fetch Position

        self.InFile.seek(self.nSeekPos[nChrom] + nFrom)            # Get Sequence

        sFetchedSeq = re.sub(re_nonchr, '', self.InFile.read(nTo - nFrom))

        print('Fetching Chromosome Sequence %s ...DONE' % sChrom)

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


def main():

    sTargetSeedType = '8mer'
    sQueue          = 'jinman.mirna3'
    bTestRun        = True
    sSegmentWindow  = 125
    sAnalysis       = 'geoES'
    sGenome         = 'mmu'


    copy_temp_core_script ()

    #Part 1: Get Segments
    qsub_Segmentation (sTargetSeedType, sGenome, sQueue, bTestRun, sSegmentWindow)

    #Part 2: Calculate Basic Stats
    #local_Basic_stats(sTargetSeedType, sSegmentWindow)

    #Part 3: Generate Dot-Bracket Notation and Calculate Min. Free Energy and LocalAU Sore
    #qsub_RNAfold_and_localAU (sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 4: Filter/Screen by Dot-Bracket Notation
    #qsub_Filter_dot_bracket (sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 5: Apply Additional Constraints
    #qsub_Apply_additional_constraints(sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 5: Compile Filtered Data,  Sort and Output Secondary Structure
    #qsub_Sort_and_evaluate_segments(sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 6: Full Compiled List, Sort and Output Secondary Structure
    #qsub_Full_compile_list (sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 7: Extend Segment Sequence for Secondary Structure
    #qsub_extended_ss (sTargetSeedType, sQueue, bTestRun, sSegmentWindow)

    #Part 8: Compare 75, 100, 125 bps list
    #qsub_compare_lists (sTargetSeedType, sQueue, bTestRun)

    #Part 9: Expression Profile from small RNA and mRNA Sequencing Data (Cell-Lines)
    #qsub_get_maturemir_RPM (sGenome, sQueue, bTestRun)

    #qsub_compile_targetmir_rpm_list (sAnalysis, sGenome, sQueue, bTestRun)

    #qsub_compile_targetmir_rpkm_list (sAnalysis, sGenome, sQueue, bTestRun)


#def END: main()


def qsub_Segmentation (sTargetSeedType, sGenome, sQueue, bTestRun, sSegmentWindow):
    sOutDir     = '%s/out/%s/%sbps/segments'     % (sBASE_DIR, sGenome, sSegmentWindow)
    sLogDir     = '%s/log/%s/%sbps/segments/%s'  % (sBASE_DIR, sGenome, sSegmentWindow, sTIME_STAMP)

    cGenome     = cFasta(dict_sGENOME_DIR[sGenome])


    for sChrID in cGenome.sChrIDList:
        for sStrand in ['+','-']:

            sScript = '%s/temp/B_core_miRDE_%s.py segmentation %s %s %s %s %s' %\
                      (sBASE_DIR, sTIME_STAMP, sOutDir, sTargetSeedType, sChrID, sStrand, sSegmentWindow)

            if bTestRun == True:
                print(sScript)
            else:
                os.makedirs(sOutDir, exist_ok=True)
                os.makedirs(sLogDir, exist_ok=True)
                os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.Segments.%s%s.%s.%sbps' %
                          (sScript, sLogDir, sQueue, sChrID, sStrand, sTargetSeedType, sSegmentWindow))
        #loop END: sStrand
    #loop END: sChrID
#def END: qsub_Segmentation


def local_Basic_stats (sTargetSeedType, sSegmentWindow):
    print('Outputing Basics Stats Table')

    dict_sSeedSeq   = load_87_conserved_mirna_seedseq()

    sInDir          = '%s/out/%sbps/additional_constraints'   % (sBASE_DIR, sSegmentWindow)
    sOutDir         = '%s/%sbps/stats' % (sBASE_DIR, sSegmentWindow)
    os.makedirs(sOutDir, exist_ok=True)

    list_sKey       = sorted(list(dict_sSeedSeq.keys()))

    sHeader1        = '\t'.join(list_sKey)
    sHeader2        = '\t'.join([dict_sSeedSeq[sKey][sTargetSeedType] for sKey in list_sKey if sKey != 'miR-203a-3p'])

    OutFile = open('%s/basic_stats_additional.txt' % sOutDir, 'w')

    OutFile.write('miRNA\t' + sHeader1 + '\n')
    OutFile.write('%s\t' % sTargetSeedType + sHeader2 + '\n')

    for sStrand in ['+', '-']:
        dict_nTotal = {}

        for sChrID in list_sCHRID_HSA:

            sChrom          = 'chr%s%s' % (sChrID, sStrand)

            list_sPrintLine = [sChrom]

            for sKey in list_sKey:

                if sKey == 'miR-203a-3p': continue

                dict_sSegments = pickle.load(open('%s/%s/%s.list' % (sInDir, sChrom, sKey), 'rb'))
                nSegments      = len(dict_sSegments)

                list_sPrintLine.append(str(nSegments))

                if sKey not in dict_nTotal:
                    dict_nTotal[sKey] = 0
                dict_nTotal[sKey] += nSegments
            #loop END: sKey

            sOutLine = '\t'.join(list_sPrintLine)

            OutFile.write(sOutLine + '\n')
        #loop END: sChrID

        sTotal  = '\t'.join([str(dict_nTotal[sKey]) for sKey in list_sKey if sKey != 'miR-203a-3p'])

        OutFile.write('Total\t' + sTotal + '\n')
    #loop END: sStrand

    OutFile.close()

    print('Outputing Basics Stats Table...DONE')
#def END: local_Basic_stats


def qsub_RNAfold_and_localAU (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir       = '%s/out/%sbps/segments'            % (sBASE_DIR, sSegmentWindow)
    sLogDir      = '%s/log/%sbps/rnafold_localAU/%s'  % (sBASE_DIR, sSegmentWindow, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)

    list_sChrIDs = ['chr%s%s' % (sChrID, sStrand) for sChrID in list_sCHRID_HSA for sStrand in ['+', '-']]

    for sChrID in list_sChrIDs:

        sOutDir  = '%s/out/%sbps/rnafold_localAU/%s' % (sBASE_DIR, sSegmentWindow, sChrID)
        os.makedirs(sOutDir, exist_ok=True)

        sScript  = '%s/temp/B_core_miRDE_%s.py RNAfold_and_localAU %s %s %s %s %s' % \
                    (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sChrID, sTargetSeedType, sSegmentWindow)

        if bTestRun == True:
            print(sScript)
        else:
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.RNAfold.%s.%s.%sbps' %
                      (sScript, sLogDir, sQueue, sChrID, sTargetSeedType, sSegmentWindow))
    #loop END: sChrID
#def END: qsub_RNAfold_and_localAU


def qsub_Filter_dot_bracket (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir       = '%s/out/%sbps/rnafold_localAU'     % (sBASE_DIR, sSegmentWindow)
    sLogDir      = '%s/log/%sbps/dotbr_filtered/%s'   % (sBASE_DIR, sSegmentWindow, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)

    list_sChrIDs = ['chr%s%s' % (sChrID, sStrand) for sChrID in list_sCHRID_HSA for sStrand in ['+', '-']]

    for sChrID in list_sChrIDs:

        sOutDir     = '%s/out/%sbps/dotbr_filtered/%s'       % (sBASE_DIR, sSegmentWindow, sChrID)
        os.makedirs(sOutDir, exist_ok=True)

        sScript     = '%s/temp/B_core_miRDE_%s.py filter_dot_bracket %s %s %s %s' % \
                        (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sChrID, sSegmentWindow)

        if bTestRun == True:
            print(sScript)
        else:
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.Filter.%s.%s.%sbps' %
                      (sScript, sLogDir, sQueue, sChrID, sTargetSeedType, sSegmentWindow))
    #loop END: sChrID
#def END: qsub_Filter_dot_bracket


def qsub_Apply_additional_constraints (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir       = '%s/out/%sbps/dotbr_filtered'              % (sBASE_DIR, sSegmentWindow)
    sLogDir      = '%s/log/%sbps/additional_constraints/%s'   % (sBASE_DIR, sSegmentWindow, sTIME_STAMP)
    os.makedirs(sLogDir, exist_ok=True)

    list_sChrIDs = ['chr%s%s' % (sChrID, sStrand) for sChrID in list_sCHRID_HSA for sStrand in ['+', '-']]

    for sChrID in list_sChrIDs:

        sOutDir     = '%s/out/%sbps/additional_constraints/%s'% (sBASE_DIR, sSegmentWindow, sChrID)
        os.makedirs(sOutDir, exist_ok=True)

        sScript     = '%s/temp/B_core_miRDE_%s.py apply_additional_constraints %s %s %s %s' % \
                        (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sChrID, sSegmentWindow)

        if bTestRun == True:
            print(sScript)
        else:
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.AddConstraints.%s.%s.%sbps' %
                      (sScript, sLogDir, sQueue, sChrID, sTargetSeedType, sSegmentWindow))
    #loop END: sChrID
#def END: qsub_Apply_additional_constraints


def qsub_Sort_and_evaluate_segments (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir  = '%s/out/%sbps/additional_constraints'   % (sBASE_DIR, sSegmentWindow)
    sOutDir = '%s/out/%sbps/sort_and_evaluate'        % (sBASE_DIR, sSegmentWindow)
    sLogDir = '%s/log/%sbps/sort_and_evaluate/%s'     % (sBASE_DIR, sSegmentWindow, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    dict_sSeedSeq = load_87_conserved_mirna_seedseq()

    list_sMirName = sorted(list(dict_sSeedSeq.keys()))


    for sMirName in list_sMirName:

        if sMirName == 'miR-203a-3p': continue

        sScript = '%s/temp/B_core_miRDE_%s.py sort_and_evaluate_segments %s %s %s %s' % \
                  (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sMirName, sSegmentWindow)

        if bTestRun == True:
            print(sScript)
        else:
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.SortEval.%s.%s.%sbps' %
                      (sScript, sLogDir, sQueue, sMirName, sTargetSeedType, sSegmentWindow))
    #loop END: sMirName
#def END: qsub_Sort_and_evaluate_segments


def qsub_Full_compile_list (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir  = '%s/out/%sbps/sort_and_evaluate'    % (sBASE_DIR, sSegmentWindow)
    sOutDir = '%s/out/%sbps/full_compile_list'    % (sBASE_DIR, sSegmentWindow)
    sLogDir = '%s/log/%sbps/full_compile_list/%s' % (sBASE_DIR, sSegmentWindow, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sScript = '%s/temp/B_core_miRDE_%s.py full_compile_list %s %s %s' % \
              (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sSegmentWindow)

    if bTestRun == True:
        print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.FullList.%s.%sbps' %
                  (sScript, sLogDir, sQueue, sTargetSeedType, sSegmentWindow))
    #if END: bTestRun
#def END: qsub_Full_compile_list


def qsub_extended_ss (sTargetSeedType, sQueue, bTestRun, sSegmentWindow):
    sInDir  = '%s/out/75bps/full_compile_list'     % sBASE_DIR
    sOutDir = '%s/out/75bps/extended_ss'           % sBASE_DIR
    sLogDir = '%s/log/75bps/extended_ss/%s'        % (sBASE_DIR, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sScript = '%s/temp/B_core_miRDE_%s.py extended_ss %s %s %s' % \
              (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir, sSegmentWindow)

    if bTestRun == True:
        print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.Extended.%s.%sbps' %
                  (sScript, sLogDir, sQueue, sTargetSeedType, sSegmentWindow))
    #if END: bTestRun
#def END: qsub_Full_compile_list


def qsub_compare_lists (sTargetSeedType, sQueue, bTestRun):

    sInDir  = '%s/out'                     % sBASE_DIR
    sOutDir = '%s/out/compare_lists'        % sBASE_DIR
    sLogDir = '%s/log/compare_lists/%s'     % (sBASE_DIR, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sScript = '%s/temp/B_core_miRDE_%s.py compare_lists %s %s' % \
              (sBASE_DIR, sTIME_STAMP, sInDir, sOutDir)

    if bTestRun == True:
        print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.CompareList.%s' %
                  (sScript, sLogDir, sQueue, sTargetSeedType))
    #if END: bTestRun
#def END: qsub_compare_lists


def qsub_get_maturemir_RPM (sAnalysis, sQueue, bTestRun):

    sInDir  = '%s/%s'                       % (sALN_DIR, sAnalysis)
    sOutDir = '%s/out/expr/sRNA'            % sBASE_DIR
    sLogDir = '%s/log/expr/%s'              % (sBASE_DIR, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    list_sFiles     = get_file_id_list(sAnalysis, 'sRNA')

    for sFileName in list_sFiles:

        sScript = '%s/temp/B_core_miRDE_%s.py get_maturemir_RPM %s %s %s %s' % \
                  (sBASE_DIR, sTIME_STAMP, sAnalysis, sFileName, sInDir, sOutDir)

        if bTestRun == True:
            print(sScript)
        else:
            os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.sRNA_RPM.%s' %
                      (sScript, sLogDir, sQueue, sFileName))
        #if END: bTestRun
    #loop END: sFileName
#def END: qsub_get_maturemir_RPM


def qsub_compile_targetmir_rpm_list (sAnalysis, sGenome, sQueue, bTestRun):

    sInDir  = '%s/out/expr/sRNA'            % sBASE_DIR
    sOutDir = '%s/out/expr/compile_list'    % sBASE_DIR
    sLogDir = '%s/log/expr/%s'              % (sBASE_DIR, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sScript = '%s/temp/B_core_miRDE_%s.py compile_targetmir_rpm_list %s %s %s' % \
              (sBASE_DIR, sTIME_STAMP, sGenome, sInDir, sOutDir)

    if bTestRun == True:
        print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.CompileListRPM' %
                  (sScript, sLogDir, sQueue))
    #if END: bTestRun
#def END: qsub_compile_targetmir_rpm_list


def qsub_compile_targetmir_rpkm_list (sAnalysis, sGenome, sQueue, bTestRun):

    sInDir  = '/extdata5/Jinman/03_CG_mrSeq/01_output/%s' % sAnalysis
    sOutDir = '%s/out/expr/compile_list'     % sBASE_DIR
    sLogDir = '%s/log/expr/%s'               % (sBASE_DIR, sTIME_STAMP)

    os.makedirs(sOutDir, exist_ok=True)
    os.makedirs(sLogDir, exist_ok=True)

    sScript = '%s/temp/B_core_miRDE_%s.py get_target_mir_rpkm_list %s %s %s' % \
              (sBASE_DIR, sTIME_STAMP, '%s_%s' % (sGenome, sAnalysis), sInDir, sOutDir)

    if bTestRun == True:
        print(sScript)
    else:
        os.system('echo "%s" | qsub -cwd -j y -o %s -q %s -N Jinman.mirDE.CompileListRPKM' %
                  (sScript, sLogDir, sQueue))
    #if END: bTestRun
#def END: qsub_compile_targetmir_rpkm_list

main()

