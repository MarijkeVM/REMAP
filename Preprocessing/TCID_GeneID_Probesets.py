
import argparse
import pandas as pd
import csv

parser = argparse.ArgumentParser(description="Extract the manufacturer's gene ID and annotated gene name from a file as well as the associated probe sets.")
parser.add_argument('TranscriptAnnotfile', metavar='TranscriptAnnotfile', type=str, nargs=1, help="The name of the transcript annotation file")
parser.add_argument('ProbesetAnnotfile', metavar='ProbesetAnnotfile', type=str, nargs=1, help="The name of the probeset annotation file")
parser.add_argument('location', metavar='location', type=str, nargs=1, help="The location of the output and indexing files")
parser.add_argument('prefix', metavar='prefix', type=str, nargs=1,help="The prefix of the output file")

args=parser.parse_args()

def main(argv):
	TrAnnotfile=argv.TranscriptAnnotfile[0]
	PrAnnotfile=argv.ProbesetAnnotfile[0]
	prefix=argv.prefix[0]
	with open(TrAnnotfile,"r") as file:
		with open(location+"/"+prefix+"_"+"TCID_GeneID.txt","w") as output:
			output.write("{0}\t{1}\n".format("TCID","GeneID"))
			for i in range(20):
				file.readline()
			for line in file:
				#print(line)
				row=line.split(",")
				TCID=row[0].strip('"')
				TCID=TCID.split(".")[0]
				GeneID=row[7].strip('"')
				GeneIDs=GeneID.split("//")
				GeneIDs=[x.strip(" ") for x in GeneIDs]
				if len(GeneIDs)>2:
					HGNC=range(1,len(GeneIDs),5)
					G=[GeneIDs[g] for g in HGNC]
					G=set(G)
					output.write("{0}\t{1}\n".format(TCID," /// ".join(G)))
				else:
					G="---"
					output.write("{0}\t{1}\n".format(TCID,G))
		
	with open(PrAnnotfile,"r") as file:
		with open(location+"/"+"ProbesetID_Name_TC.txt","a") as outfile:
			for line in file:	
				if line[1]!="J" and line[1]!="P":
					continue			
				line=line.split(",")
				PSR=line[0].strip('\"')
				PSR=PSR.split(".")[0]
				TC=line[6].strip('\"')
				TC=TC.split("///")
				TCs=list(set(TC))
				Tc=TC[0]
				Tc=Tc.split(".")[0]
				row1=[PSR,Tc]
				outfile.write("{0}\t{1}\n".format(row1[0],row1[1]))

	File=pd.read_table(location+"/"+"ProbesetID_Name_TC.txt",header=None)
	File=File.sort_values(by=[1,0])
	File.to_csv(location+"/"+"ProbesetID_Name_TC.txt",sep="\t",header=False,index=False)	
	
	DF2=pd.read_table(location+"/"+"ProbesetID_Name_TC.txt",sep="\t",header=None)

	with open(location+"/"+prefix+"_"+"TCID_GeneID.txt","r") as input1:
		header=input1.readline()
		with open(location+"/"+prefix+"_"+"TCID_GeneID_Probesets.txt","w") as out2:
			for line1 in input1:
				line=line1.split("\t")
				TC=line[0]
				Gene=line[1]
				PSRs=DF2.ix[DF2.ix[:,1]==TC,0]
				DF2=DF2.drop(DF2.index[PSRs.index])
				R=[TC]+[Gene]+["|".join(PSRs.tolist())]
				R[1]=R[1].rstrip("\n")
				out2.write("{0}\t{1}\t{2}\n".format(R[0],R[1],R[2]))
		
	#DF1=pd.read_table(prefix+"TCID_GeneID.txt",sep="\t",header=None)
	#DF2=pd.read_table("HTA-2_0_ProbesetID_Name_TC.txt",sep="\t",header=None)	
	
	
	# for index, row in DF1.iterrows():
		#print(row)
		# PSRs=DF2.ix[DF2.ix[:,1]==row[0],0].tolist()
		# TCinfo=row.values.tolist()
		# if len(PSRs)>0:
			# R=TCinfo+["|".join(PSRs)]
		# else:
			# R=TCinfo+PSRs
		# with open(prefix+"TCID_GeneID_ProbeSets.txt","a") as out2:
			# writer = csv.writer(out2, delimiter="\t")
			# writer.writerow(R)
		# out2.close()	
		
						
if __name__ == "__main__":
	main(args)