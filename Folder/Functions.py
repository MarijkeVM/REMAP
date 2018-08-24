# Code to merge ensembl gene information with annotation of probes and junctions to get an idea of the gene structure

import pandas as pd
import re
from collections import Counter
from itertools import permutations
import csv
import string 
from cogent.db.ensembl import Genome, HostAccount
from difflib import SequenceMatcher   

#account = HostAccount('ensembldb.ensembl.org', 'anonymous', '', port=5306)
account = HostAccount('127.0.0.1', 'root', 'ensembl', port=3306)
Release = 89
HumanDB = Genome(Species='human', Release=Release, account=account)

class SettingTCID_GeneIDError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

def GeneStructure(TCID,GeneID,Probesets=None,MappingFile="HTA_2_0_Probeset_SequenceIndices.txt",SequenceFile="HTA_2_0_Probeset_Sequences.txt",location="Output",prefix=""):
    
    if GeneID=="---":
        out=[TCID, GeneID, "No Annotated Gene"]
        return out

    Map=pd.read_table(MappingFile,header=0)
    GeneIDs=GeneID.split(" /// ")
    values=GeneIDs
    #print values
    #print len(Probesets)
    #TCID = 'TC1500264'
    #values = ['CAPN3', 'GANC']

    GeneDict = {'TCID_GeneID': [], 'Sequence': [], 'CDS': [], 'Exons': [], 'Transcript': [], 'Chromosome': [],
                'Start': [],
                'End': [], 'Strand': [], 'ENSG': [], 'GeneID': []}
    ExonDict = {'ENSE': [], 'Sequence': [], 'ENSG': [], 'GeneID': []}
    IntronDict = {'ENSE': [], 'Sequence': [], 'ENSG': [], 'GeneID': []}
    for v in values:
        v=v.strip()
        genes = HumanDB.getGenesMatching(Symbol=v.upper())
        #print v
        for gene in genes:
            if gene.Symbol.upper() == v.upper():
                break
        try:
            notfound=gene.Symbol
        except UnboundLocalError:
            GeneTranscripts="---"
            Annotated=[] 
            NotAnnotated=[]
            JUCAnnotated=[]
            JUCNotAnnotated=[]
            UniqueAnnotatedGenes =[]
            return GeneTranscripts, Annotated, NotAnnotated, JUCAnnotated,JUCNotAnnotated, UniqueAnnotatedGenes  
        for tr in gene.Transcripts:
            #print(tr)
            GeneDict['TCID_GeneID'].append(TCID + '_' + v)
            GeneDict['Sequence'].append(str(tr.Seq))
            GeneDict['CDS'].append(str(tr.Cds))
            E = [e.StableId for e in tr.Exons]
            GeneDict['Exons'].append('|'.join(E))
            GeneDict['Transcript'].append(tr.StableId)
            Loc = tr.Location
            Locs = str(Loc).split(":")
            GeneDict['Chromosome'].append(Locs[2])
            StEn = Locs[3].split('-')
            GeneDict['Start'].append(StEn[0])
            GeneDict['End'].append(StEn[1])
            GeneDict['Strand'].append(Locs[-1])
            GeneDict['ENSG'].append(gene.StableId)
            GeneDict['GeneID'].append(v)

            for e1 in tr.Exons:
                ExonDict['ENSE'].append(e1.StableId)
                ExonDict['Sequence'].append(str(e1.Seq))
                ExonDict['ENSG'].append(gene.StableId)
                ExonDict['GeneID'].append(v)
                
            if tr.Introns is not None:
                for i1 in tr.Introns:
                    IntronDict['ENSE'].append('Intron_' + str(i1.Rank) + "_" + gene.StableId)
                    IntronDict['Sequence'].append(str(i1.Seq))
                    IntronDict['ENSG'].append(gene.StableId)
                    IntronDict['GeneID'].append(v)

    GeneTranscripts = pd.DataFrame(data=GeneDict)
    #print GeneTranscripts
    GeneTranscripts = GeneTranscripts.reindex(
        ['TCID_GeneID', 'Sequence', 'Exons', 'Transcript', 'Chromosome', 'Start', 'End', 'Strand', 'ENSG', 'GeneID','CDS'],
        axis=1)
    ExonInfo1 = pd.DataFrame(data=ExonDict)
    ExonInfo1 = ExonInfo1.reindex(['ENSE', 'Sequence', 'ENSG', 'GeneID'], axis=1)
    ExonInfo1 = ExonInfo1.drop_duplicates()
    IntronInfo1 = pd.DataFrame(data=IntronDict)
    IntronInfo1 = IntronInfo1.reindex(['ENSE', 'Sequence', 'ENSG', 'GeneID'], axis=1)
    IntronInfo1 = IntronInfo1.drop_duplicates()

    ExonInfo = pd.concat([ExonInfo1, IntronInfo1])
    ExonInfo.index = range(len(ExonInfo.index))
    #ExonInfo.to_csv("Exoninfo.csv",sep="\t")
    PSR = [p for p in Probesets if 'PSR' in p]
    Annotated = []
    NotAnnotated = []
    for p in PSR:
        # p = p.split("_")[0]
        print p
        Rows = pd.DataFrame(Map.loc[Map['ProbesetID'] == p])
        # print Rows
        Rows.index = range(1)
        columns = pd.read_table(SequenceFile, header=0,nrows=1).columns
        ProbeInfo = pd.read_table(SequenceFile, header=None, skiprows=Rows.ix[0, 1] - 1, nrows=Rows.ix[0, 2],
                                  usecols=["probeset_id", "probe_id", "probe_sequence"],names=columns)

        ProbeInfo = ProbeInfo.sort_values(by="probe_id")
		#print(ProbeInfo)
        ProbeInfo.index = range(len(ProbeInfo.index))
        for s in range(len(ProbeInfo.index)):
            if s == 0:
                Seq = ProbeInfo.ix[s, "probe_sequence"]
                N = str(ProbeInfo.ix[s, "probe_id"])
            else:
                # ExonPosEndPrev = ProbeInfo.ix[s - 1, 4] + len(ProbeInfo.ix[s - 1, 9])
                # ExonPosEnd = ProbeInfo.ix[s, 4] + len(ProbeInfo.ix[s, 9])
                # Rest = ExonPosEnd - ExonPosEndPrev
                # if Rest < 25:
                    # Seq += ProbeInfo.ix[s, 9][len(ProbeInfo.ix[s, 9]) - Rest:]
                    # N += "_" + str(ProbeInfo.ix[s, 5])
                # else:
                Seq += "..." + ProbeInfo.ix[s, "probe_sequence"]
                N += "..." + str(ProbeInfo.ix[s, "probe_id"])
        Pieces = Seq.split("...")
        basecomplement = string.maketrans('ACGT', 'TGCA')
        Temp = [v.translate(basecomplement) for v in Pieces]
        Pieces = Temp

        E = pd.DataFrame(index=Pieces, columns=['Match','ENSG','GeneID'])
        for pi in Pieces:
            r = re.compile(pi)
            Match = ""
            ENSGExon= ""
            GNameExon=""
            for exon in range(len(ExonInfo.index)):
                Find = r.search(ExonInfo.ix[exon, 1])
                if Find != None:
                    if len(Match) == 0:
                        Match += ExonInfo.ix[exon, 0]
                        ENSGExon += ExonInfo.ix[exon, 2]
                        GNameExon += ExonInfo.ix[exon, 3]
                    else:
                        Match += "|" + ExonInfo.ix[exon, 0]
                        ENSGExon += "|"+ExonInfo.ix[exon, 2]
                        GNameExon += "|"+ExonInfo.ix[exon, 3]
            E.ix[pi,0] = Match
            E.ix[pi,1] = ENSGExon
            E.ix[pi,2] = GNameExon

        if all(E.ix[pi, 0] == "" for pi in Pieces):
            NotAnnotated.append([TCID, p, "-", "-", "-"])
            continue

        MatchedExons = pd.DataFrame({'Match': list(E.ix[:, 0]), "Probes": N.split("..."), "ENSG": list(E.ix[:, 1]),"GeneID": list(E.ix[:, 2])})
        MatchedExons = MatchedExons[['Match','Probes','ENSG','GeneID']]
        #print MatchedExons
        #MatEx = set(MatchedExons.ix[:, 0])
        # if len(MatEx) == 1:  # all probes match the same exons
            # if len(MatchedExons.index) == len(Pieces):
                # Sep = MatchedExons.ix[0, 0]
                # L = Sep.split("|")
                #print L
                # Genes = ExonInfo.loc[ExonInfo.ix[:, 0].isin(L), ['ENSG', 'GeneID']]
                #Set1 = list(Genes.ix[:, 0])
                #print Set1
                #SetI = "_".join(Set1)
                #Set2 = list(Genes.ix[:, 1])
                #print Set2
                #SetII = "_".join(Set2)
                # GNames=list(set(Genes.ix[:, 1]))
                # GNamesS="_".join(GNames)
                # Annotated.append([TCID + '_' + GNamesS, p, N, Seq, MatchedExons.ix[0, 0], MatchedExons.ix[0,2], MatchedExons.ix[0,3]])
            # else:
                # Sep = MatchedExons.ix[0, 0]
                # L = Sep.split("|")
                # Genes = ExonInfo.loc[ExonInfo.ix[:, 0].isin(L), ['ENSG', 'GeneID']]
                #Set1 = list(Genes.ix[:, 0])
                #print Set1
                #SetI = "_".join(Set1)
                #Set2 = list(Genes.ix[:, 1])
                #print Set2
                #SetII = "_".join(Set2)
                # GNames=list(set(Genes.ix[:, 1]))
                # GNamesS="_".join(GNames)
                # PartiallyAnnotated.append([TCID + '_' + GNamesS, p, N, Seq, MatchedExons.ix[0, 0],MatchedExons.ix[0,2], MatchedExons.ix[0,3]])
        #else:
				
        All = "|".join(MatchedExons.ix[:, 0])
        AllENSG="|".join(MatchedExons.ix[:, 2])
        AllG="|".join(MatchedExons.ix[:, 3])
        Sep = All.split("|")
        Counts = Counter(Sep)
			
        AnnotatedByAllProbes_temp = [x for x in Counts.keys() if Counts[x]%len(Pieces)==0]
        AnnotatedBySubsetProbes = [x for x in Counts.keys() if Counts[x] < len(Pieces)]
           
        AnnotatedByAllProbes=[]
        for x in AnnotatedByAllProbes_temp:
            C=0
            for r in range(len(MatchedExons.index)):
               if x in MatchedExons.ix[r,0]:
                   C+=1
            if C==len(Pieces):
                AnnotatedByAllProbes.append(x)				
            else:
                AnnotatedBySubsetProbes.append(x)

        ENSGSep=AllENSG.split("|")		
        GSep=AllG.split("|")
		
        ETemp=pd.DataFrame()
        ETemp['All']=Sep
        ETemp['ENSG']=ENSGSep
        ETemp['G']=GSep
        ETemp=ETemp.drop_duplicates(keep='first') 
			
            #ExonsI = "|".join(AnnotatedByAllProbes)
            #ExonsII = "|*".join(AnnotatedBySubsetProbes)
            # if len(ExonsI)>0 and len(ExonsII)>0:
                # Exons = ExonsI + "|*" + ExonsII
            # elif len(ExonsI)>0 and len(ExonsII)==0:
                # Exons = ExonsI
            # elif len(ExonsI)==0 and len(ExonsII)>0:
                # Exons="*"+ExonsII
            # else:
                # Exons="-"			

        L = AnnotatedByAllProbes+AnnotatedBySubsetProbes
        ExonsI=""
        SetIa=""
        SetIIa=""
        for l in AnnotatedByAllProbes:
            settemp=ETemp.loc[ETemp.ix[:,0]==l,:]
            settemp1=settemp['All'].tolist()
            settemp2=[x.split("_")[0] for x in settemp1]
            ExonsI=ExonsI+"|".join(settemp2)+"|"
            SetIa=SetIa+"|".join(settemp.ix[:, 1])+"|" 
            SetIIa=SetIIa+"|".join(settemp.ix[:, 2])+"|" 	
        ExonsI=ExonsI.rstrip("|")
        SetIa=SetIa.rstrip("|")
        SetIIa=SetIIa.rstrip("|")			
        ExonsII=""
        SetIb=""
        SetIIb=""		
        for l in AnnotatedBySubsetProbes:
            settemp=ETemp.loc[ETemp.ix[:,0]==l,:]
            settemp1=settemp['All'].tolist()
            settemp2=[x.split("_")[0] for x in settemp1]
            ExonsII=ExonsII+"|*".join(settemp2)+"|*"
            SetIb=SetIb+"|*".join(settemp.ix[:, 1])+"|*" 
            SetIIb=SetIIb+"|*".join(settemp.ix[:, 2])+"|*" 
        ExonsII=ExonsII.rstrip("|*")
        SetIb=SetIb.rstrip("|*")
        SetIIb=SetIIb.rstrip("|*")
        if len(ExonsI)>0 and len(ExonsII)>0:
            Exons = ExonsI + "|*" + ExonsII
            SetI=SetIa+"|*"+SetIb
            SetII=SetIIa+"|*"+SetIIb			
        elif len(ExonsI)>0 and len(ExonsII)==0:
            Exons = ExonsI
            SetI=SetIa
            SetII=SetIIa
        elif len(ExonsI)==0 and len(ExonsII)>0:
            Exons="*"+ExonsII
            SetI="*"+SetIb
            SetII="*"+SetIIb			
        else:
            Exons="-"		
			
            #print L
            #print ExonInfo.loc[ExonInfo.ix[:, 0].isin(L), ['ENSG', 'GeneID']]			
        Genes = ExonInfo.loc[ExonInfo.ix[:, 0].isin(L), ['ENSG', 'GeneID']]
            #print Genes
            #Set1 = list(Genes.ix[:, 0])
            #SetI = "_".join(Set1)
            #Set2 = list(Genes.ix[:, 1])
            #SetII = "_".join(Set2)
        GNames=list(set(Genes.ix[:, 1]))
        GNamesS="_".join(GNames)
        Annotated.append([TCID + '_' + GNamesS, p, Exons, SetI, SetII])

    # The remaining probe sets also belong to a gene: this can be determined by the order of the probe sets
    AllocatedProbes = [x[1] for x in Annotated]
    All = Annotated
    print NotAnnotated
    NotAnnotatedtemp=NotAnnotated
    for e in range(len(NotAnnotatedtemp)):
        #print e
        PSRtemp = NotAnnotatedtemp[e][1]
        Index = PSR.index(PSRtemp)
        IndexPrev = Index - 1
        TCPrev = None
        TCNext = None

        if IndexPrev >= 0:
            while (not PSR[IndexPrev] in AllocatedProbes) & (IndexPrev >= 0):
                IndexPrev -= 1
				
        IndexNext = Index + 1
        if IndexNext <= (len(PSR) - 1):
            while (not PSR[IndexNext] in AllocatedProbes) & (IndexNext <= (len(PSR) - 1)):
                IndexNext += 1
                if (IndexNext > (len(PSR) - 1)):
                    break
				   
        if IndexPrev >= 0:
            TCPrev = [x[0] for x in All if x[1] == PSR[IndexPrev]][0]
        if IndexNext <= (len(PSR) - 1):
            TCNext = [x[0] for x in All if x[1] == PSR[IndexNext]][0]

        #if (IndexPrev < 0) & (TCNext is not None):
        #    TCPrev = TCNext
        #elif (TCPrev is not None) & (IndexNext > (len(PSR) - 1)):
        #    TCNext = TCPrev
        #if len(NotAnnotated) == len(PSR):  # none of the psr is annotated
        #    TCPrev = TCID + "_" + "_".join(GeneIDs)
        #    TCNext = TCPrev

        if TCNext == TCPrev:
            NotAnnotated[e][0] = TCPrev
        # else:
            # t = NotAnnotated[e]
            # t = t[1:]
            # t1 = [TCPrev] + t
            # t2 = [TCNext] + t
            # NotAnnotated[e] = t1
            # NotAnnotated.insert(e + 1, t2)

    JUC = [j for j in Probesets if "JUC" in j]
    JUCAnnotated = []
    JUCNotAnnotated = [] 
    for j in JUC:
        # j = j.split("_")[0]
        print j
        Rows = pd.DataFrame(Map.loc[Map['ProbesetID'] == j])
        Rows.index = range(1)
        columns = pd.read_table(SequenceFile, header=0,nrows=1).columns
        ProbeInfo = pd.read_table(SequenceFile, header=None, skiprows=Rows.ix[0, 1] - 1, nrows=Rows.ix[0, 2],
                                  usecols=["probeset_id", "probe_id", "probe_sequence"],names=columns)

        ProbeInfo = ProbeInfo.sort_values(by="probe_id")
        ProbeInfo.index = range(len(ProbeInfo.index))
        StartPlace=list()
        EndPlace=list()
        for s in range(1,len(ProbeInfo.index)):
            if s<=1:
                #Seq = ProbeInfo.ix[s, 9]
                #N = str(ProbeInfo.ix[s, 4])
                str1=ProbeInfo.ix[0, "probe_sequence"]
                str2=ProbeInfo.ix[1, "probe_sequence"]
                match = SequenceMatcher(None, str1, str2).find_longest_match(0, len(str1), 0, len(str2))
                Commonstr=str1[match.a: match.a + match.size]
            else:
                #ExonPosEndPrev = ProbeInfo.ix[s - 1, 4] + len(ProbeInfo.ix[s - 1, 9])
                # ExonPosEnd = ProbeInfo.ix[s, 4] + len(ProbeInfo.ix[s, 9])
                # Rest = ExonPosEnd - ExonPosEndPrev
                # if Rest < 25:
                    # Seq += ProbeInfo.ix[s, 9][len(ProbeInfo.ix[s, 9]) - Rest:]
                    # N += "_" + str(ProbeInfo.ix[s, 5])
                # else:
                #Seq += "..." + ProbeInfo.ix[s, 9]
                #N += "..." + str(ProbeInfo.ix[s, 4])
                str3=ProbeInfo.ix[s, "probe_sequence"]
                match = SequenceMatcher(None, str3, Commonstr).find_longest_match(0, len(str3), 0, len(Commonstr))
                Commonstr=str3[match.a: match.a + match.size]
        for s in range(len(ProbeInfo.index)):       
                str1=ProbeInfo.ix[s, "probe_sequence"]
                StartPlace.append(str1.index(Commonstr))
                EndPlace.append(str1.index(Commonstr)+len(Commonstr) - 1)
		
        First=StartPlace.index(max(StartPlace))
        Last=EndPlace.index(min(EndPlace))
		
        Seq=ProbeInfo.ix[First, "probe_sequence"][0:max(StartPlace)]+Commonstr+ProbeInfo.ix[Last, "probe_sequence"][min(EndPlace)+1:]

		
        Piece = Seq
        basecomplement = string.maketrans('ACGT', 'TGCA')
        Temp = Piece.translate(basecomplement)
        Piece = Temp

        r = re.compile(Piece)
        FoundM = []
        FoundS = []
        GeneIM = []
        GeneIIM = []
        GeneIS = []
        GeneIIS = []		
        for transcript in range(len(GeneTranscripts.index)):
            Find = r.search(GeneTranscripts.ix[transcript, 10])
            if Find != None:
                Exons = GeneTranscripts.ix[transcript, 2].split("|")
                GeneF = GeneTranscripts.ix[transcript, 9]
                GeneENSG = GeneTranscripts.ix[transcript, 8]
                C_temp = permutations(Exons, 2)
                C = [c_temp for c_temp in C_temp]
                for s in C:
                    Index1 = ExonInfo.loc[ExonInfo.ix[:, 0] == s[0], :].index[0]
                    Index2 = ExonInfo.loc[ExonInfo.ix[:, 0] == s[1], :].index[0]
                    PastedSeq = ExonInfo.ix[Index1, 1] + ExonInfo.ix[Index2, 1]
                    Match = r.search(PastedSeq)
                    Comb = s[0]+"_"+GeneF + "|" + s[1]+"_"+GeneF
                    if Match != None and Comb not in FoundM:
                        FoundM.append(Comb)
                        GeneIM.append(GeneF+"|"+GeneF)
                        GeneIIM.append(GeneENSG+"|"+GeneENSG)
            else:
                Exons = GeneTranscripts.ix[transcript, 2].split("|")
                GeneF = GeneTranscripts.ix[transcript, 9]
                GeneENSG = GeneTranscripts.ix[transcript, 8]
                for s in Exons:
                    Index = ExonInfo.loc[ExonInfo.ix[:, 0] == s, :].index[0]
                    SeqE = ExonInfo.ix[Index, 1]
                    Piece1 = Piece[0:(int(len(Piece) / 2))]
                    Piece2 = Piece[(int(len(Piece) / 2)):]
                    r1 = re.compile(Piece1)
                    r2 = re.compile(Piece2)
                    M1 = r1.search(SeqE)
                    M2 = r2.search(SeqE)
                    if ((M1 != None) | (M2 != None)) and (s+"_"+GeneF not in FoundS):
                        FoundS.append(s+"_"+GeneF)
                        GeneIS.append(GeneF)
                        GeneIIS.append(GeneENSG)

        GNamesJ1="|".join(GeneIM)+"|"+"|".join(GeneIS)
        GNamesJ2=GNamesJ1.split("|")
        GNamesJ2=[x for x in GNamesJ2 if len(x)>0]
        #print GNamesJ2
        #print  "_".join(set(GNamesJ2))
        if len(FoundM) > 0 and len(FoundS)==0:
            FoundMtemp1=[re.split("\|",x) for x in FoundM]
            print FoundMtemp1 
            FoundMtemp2=sum(FoundMtemp1,[])
            FoundMtemp3=[x.split("_")[0] for x in FoundMtemp2]
            FoundMkeep=[x+"|"+y for x,y in zip(FoundMtemp3[::2],FoundMtemp3[1::2])]
            JUCAnnotated.append(
                [TCID + "_" + "_".join(set(GNamesJ2)), j, "//".join(FoundMkeep), "//".join(GeneIIM), "//".join(GeneIM)])
        elif len(FoundS) > 0 and len(FoundM)==0:
            FoundSkeep=[x.split("_")[0] for x in FoundS]
            JUCAnnotated.append([TCID + "_" + "_".join(set(GNamesJ2)), j, "*"+"//*".join(FoundSkeep), "*"+"//*".join(GeneIIS), "*"+"//*".join(GeneIS)])
        elif len(FoundS) > 0 and len(FoundM)>0:	
            FoundSkeep=[x.split("_")[0] for x in FoundS]
            FoundMtemp1=[re.split("\|",x) for x in FoundM]
            #print FoundMtemp1 
            FoundMtemp2=sum(FoundMtemp1,[])
            FoundMtemp3=[x.split("_")[0] for x in FoundMtemp2]
            FoundMkeep=[x+"|"+y for x,y in zip(FoundMtemp3[::2],FoundMtemp3[1::2])]			
            JUCAnnotated.append(
                [TCID + "_" + "_".join(set(GNamesJ2)), j, "//".join(FoundMkeep)+"//*"+"//*".join(FoundSkeep), "//".join(GeneIIM)+"//*"+"//*".join(GeneIIS), "//".join(GeneIM)+"//*"+"//*".join(GeneIS)])
        elif len(FoundS) == 0 and len(FoundM)==0:
            JUCNotAnnotated.append([TCID, j, "-", "-", "-"])

    JuncsMultipleGenes = []
    for ju in range(len(JUCAnnotated)):
        J = JUCAnnotated[ju]
        genes = J[0].split("_")
        if len(genes) > 3:
            JuncsMultipleGenes.append(J)

    with open(location+'/'+prefix+'_JuncsMultipleGenes.txt', 'a') as output:
        writer = csv.writer(output, delimiter="\t")
        if len(JuncsMultipleGenes) == 1:
            writer.writerows(JuncsMultipleGenes[0])
        else:
            writer.writerows(JuncsMultipleGenes)


    def getKey(item):
        return item[0]

    Annotated.sort(key=getKey)
    NotAnnotated.sort(key=getKey)
    JUCAnnotated.sort(key=getKey)
    JUCNotAnnotated.sort(key=getKey)   
        
    return GeneTranscripts,Annotated,NotAnnotated,JUCAnnotated,JUCNotAnnotated            