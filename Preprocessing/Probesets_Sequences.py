#!/usr/bin/python

import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Merging the manufacturer's probe set ids, the sequence, the clf and pgf information")
parser.add_argument('clf', metavar='clf', type=str, nargs=1, help="The name of the clf file")
parser.add_argument('pgf', metavar='pgf', type=str, nargs=1, help="The name of the pgf file")
parser.add_argument('skipnrows', metavar='skipnrows', type=int, nargs=1, help="The number of lines starting with ""#"" in the .clf file and thus to skip. Default is 11.",default=11)
parser.add_argument('location', metavar='location', type=str, nargs=1, help="The location of the output and indexing files")
parser.add_argument('Outfile1', metavar='Outfile1', type=str, nargs=1, help="The name of the output file with the sequence information")
parser.add_argument('Outfile2', metavar='Outfile2', type=str, nargs=1,help="The name of the output file for the indexing")

args=parser.parse_args()


def main(argv):
    clf=argv.clf[0]
    pgf=argv.pgf[0]
    out1=argv.Outfile1[0]
    out2=argv.Outfile2[0]
    skips=argv.skipnrows[0]
    probe_x_y=pd.read_table(clf,skiprows=skips,sep="\t",header=None,index_col=False)
    probe_x_y.columns = ['probe_id', 'x','y']   
    with open(pgf, 'r') as pgffile:
        with open(location+"/"+out1,'w') as out1file:
            PSR=''
            Pos=''
            H=1
            Heads=[]
            for line in pgffile:
                if "#%header" in line:
                    linesplit=line.split('\t')
                    linesplit=[x.strip("\n") for x in linesplit]
                    h1=linesplit[0][10:]
                    if len(h1)>0:
                        Heads.append(h1)
                    Heads=Heads+linesplit[1:]
                elif line[0]=="#":
                    next
                else:
                    if H==1:
                        H=0
                        F=''
                        Heads=Heads+["x","y"]
                        for x in range(len(Heads)):
                            if x==len(Heads)-1:
                                F=F+'{'+str(x)+'}\n'
                            else:	
                                F=F+'{'+str(x)+'}\t'
                        print(Heads)
                        out1file.write(F.format(*Heads))
                    linesplit=line.split('\t')                   
                    if linesplit[0]=='':
                        if linesplit[1]=='':
                            #print(linesplit[2])
                            x_y=probe_x_y.loc[probe_x_y.ix[:,0]==int(linesplit[2]),["x","y"]].values[0]
                            Info='\t'.join(linesplit[2:])
                            Info=Info[:-1]
                            out1file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(PSR,Pos,Info,x_y[0],x_y[1]))
                        else:
                            Pos='\t'.join(linesplit[1:])
                            Pos=Pos[:-1]
                    else:
                        PSR='\t'.join(linesplit)
                        PSR=PSR[:-1]  
        out1file.close()
    with open(location+"/"+out1,'r') as out1file:
        with open(location+"/"+out2,'w') as out2file:
            out2file.write('{0}\t{1}\t{2}\n'.format('ProbesetID','Begin','NRows'))
            header=out1file.readline()
            line1=out1file.readline()
            row=line1.split()
            print(row)
            ID=row[2].split(".")[0]    
            prevID=ID                
            count=1
            begin=2
            for line in out1file:
                row=line.split()
                ID=row[2].split(".")[0]    
                if ID==prevID :
                    count +=1
                    continue
                else:
                    out2file.write('{0}\t{1}\t{2}\n'.format(prevID,begin,count))
                    begin+=count
                    count=1
                    prevID=ID
            out2file.write('{0}\t{1}\t{2}\n'.format(ID,begin,count))   
        out2file.close()
    out1file.close()
                    
if __name__ == "__main__":
    main(args)
    
       
    