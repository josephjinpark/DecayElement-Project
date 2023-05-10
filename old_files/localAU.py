def site_8m(sMiRNA) :
    return [site_7m8(sMiRNA)[0] + "A"]

def site_7m8(sMiRNA) :
    return [sMiRNA[-8:][:-1]]

def site_7a1(sMiRNA) :
    return [site_6m(sMiRNA)[0] + "A"]
def site_6m(sMiRNA) :
    return [sMiRNA[-7:][:-1]]
def site_6a1(sMiRNA):
    return [sMiRNA[-6:][:-1] + "A"]
def site_5m(sMiRNA) :
    return [sMiRNA[-6:][:-1]]
def site_o7m(sMiRNA) :
    return [sMiRNA[-9:][:-2]]
def site_o6m(sMiRNA) :
    return [sMiRNA[-8:][:-2]]

	
def SetAUweight(nSitelen):

    listWeight = []
    nStart = 2
    for n in range(0, 30) :
        listWeight.append( float(1.0)/float(nStart+n)  )
    listWeight.reverse()

    for n in range(0, nSitelen) :
        listWeight.append( 0 )

    nStart = 2
    #listWeight.append( float(1.0)/float(nStart) )
    for n in range(0, 30) :
        listWeight.append( float(1.0)/float(nStart+n) )


		
def pos_based_localAU(postart,posend,useq):

    fTotalAUscore = 0.0
    sitelen = posend-postart+1
    AUweightlist = SetAUweight(sitelen)
    if postart < 30:
        nSeqOffset = 30 - postart
    else:
        nSeqOffset = 0
    listAUWeight = AUweightlist[nSeqOffset:]
    nSeqFrom = postart - 30 + nSeqOffset
    nSeqTo = nSeqFrom + len(listAUWeight)
    sSubSeq = useq[nSeqFrom:nSeqTo]
    fAUScore = 0.0
    fMaxval = 0.0
    for nSeqPos in range(0, len(sSubSeq)):
        if sSubSeq[nSeqPos] == 'A' or sSubSeq[nSeqPos] == 'U' or sSubSeq[nSeqPos] == 'T' :
            fAUScore += listAUWeight[nSeqPos]
        fMaxval += listAUWeight[nSeqPos]
    fAUScore /= fMaxval
    return fAUScore
