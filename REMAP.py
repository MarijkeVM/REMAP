#!/usr/bin/env python

from Functions import Functions
import csv
import string
import argparse

parser = argparse.ArgumentParser(description="Execution of the REMAP function in which the annotation of probe sets is checked against the Ensembl database")
parser.add_argument('TC_GeneID_Probesets', metavar='input', type=str, nargs=1, help="The name of the "'TC_GeneID_Probesets'" file. The structure of the file is described in README.txt")
parser.add_argument('Probeset_SequenceIndices', metavar='input', type=str, nargs=1, help="The name of the "'Probeset_SequenceIndices'" file. The structure of the file is described in README.txt")
parser.add_argument('Probeset_Sequences', metavar='input', type=str, nargs=1, help="The name of the "'Probeset_Sequence'" file. The structure of the file is described in README.txt")
parser.add_argument('location', metavar='location', type=str, nargs=1, help="The location of the output folder")
parser.add_argument('prefix', metavar='input', type=str, nargs=1, help="The prefix for the output files")


args=parser.parse_args()
File=args.TC_GeneID_Probesets[0]
SeqIndFile=args.Probeset_SequenceIndices[0]
SeqFile=args.Probeset_Sequences[0]
location=argv.location[0]
prefix=args.prefix[0]

NotAnnotatedTC=[]
NotAnnotatedGeneEnsembl=[]

with open(File) as Input:
    for i, line in enumerate(Input):
        if i%500==0:
            print(str(i)+" Completed")
        TCLine=line

        TCLine=TCLine.rstrip()
        Parameters=TCLine.split("\t")
        TCID=Parameters[0]
        GeneID=Parameters[1]
        Probesets=[x for x in Parameters[2].split("|")]

        with open(location+'/'+prefix+'_IDS.txt','a') as output:
            output.write(TCID+"\n")  
    
   
        if(GeneID!='---'):
            GeneInfo, Annotated, NotAnnotated, JUCAnnotated,JUCNotAnnotated=Functions.GeneStructure(TCID=TCID,GeneID=GeneID,Probesets=Probesets,MappingFile=SeqIndFile,SequenceFile=SeqFile,location=location,prefix=prefix)
            if type(GeneInfo)!=str: 
                with open(location+'/'+prefix+'_GeneInfo.csv','a') as output:
                    GeneInfo.to_csv(output,header=False,index=False)  
                with open(location+'/'+prefix+'_PSRAnnotated.txt','a') as output1:
                    writer = csv.writer(output1,delimiter="\t")
                    writer.writerows(Annotated)
                with open(location+'/'+prefix+'_PSRNotAnnotated.txt','a') as output3:
                    writer = csv.writer(output3,delimiter="\t")
                    writer.writerows(NotAnnotated)      
                with open(location+'/'+prefix+'_JUCAnnotated.txt','a') as output4:
                    writer = csv.writer(output4,delimiter="\t")
                    writer.writerows(JUCAnnotated)
                with open(location+'/'+prefix+'_JUCNotAnnotated.txt','a') as output6:
                    writer = csv.writer(output6,delimiter="\t")
                    writer.writerows(JUCNotAnnotated) 

            else:
                NotAnnotatedGeneEnsembl.append(TCID)
                GeneInfo=[]
                Annotated=[] 
                NotAnnotated=[]
                JUCAnnotated=[]
                JUCNotAnnotated=[]
        else:
            NotAnnotatedTC.append(TCID)
            GeneInfo=[]
            Annotated=[] 
            NotAnnotated=[]
            JUCAnnotated=[]
            JUCNotAnnotated=[]
        
        with open(location+'/'+prefix+"_GeneAllInfo.txt","a+") as outputall:
            outputall.seek(0)
            header=outputall.readline()
            if header=='':
                outputall.write("TCID_GeneID\tProbeset\tAnnotationExon\tAnnotationG_a\tAnnotationG_b\r\n")
            writer = csv.writer(outputall,delimiter="\t")
            Lines=Annotated+NotAnnotated+JUCAnnotated+JUCNotAnnotated
            def getKey(item):
                return item[0]
            Lines.sort(key=getKey)    
            writer.writerows(Lines)  
            
            Genes=list(set([x[0] for x in Annotated+NotAnnotated+JUCAnnotated+JUCNotAnnotated]))
            # if len(Genes)==1:
                # PSR=[x for x in Probesets if "PSR" in x]
                # Found=[x[1] for x in Annotated]
                # numberfound=[Genes[0],len(Found),len(PSR)-len(Found),len(PSR)]
                # with open(location+'/'+prefix+"_Event.txt","a") as eventfile:
                    # writer = csv.writer(eventfile,delimiter="\t")        
                    # writer.writerow(numberfound)
                    
        with open(location+'/'+prefix+'_IDS_Completed.txt','a') as output:
            output.write(TCID+" Completed\n")  
    
if(len(NotAnnotatedGeneEnsembl)>0):    
    with open(location+'/'+prefix+"_NotAnnotatedGeneEnsembl.txt","w") as outputNotFoundGene:
        writer = csv.writer(outputNotFoundGene,delimiter="\n")
        writer.writerows([NotAnnotatedGeneEnsembl])      
    
if(len(NotAnnotatedTC)>0):    
    with open(location+'/'+prefix+"_NotAnnotatedTC.txt","w") as outputNotFoundTC:
        writer = csv.writer(outputNotFoundTC,delimiter="\n")
        writer.writerows([NotAnnotatedTC])     
                        

    
      
   
