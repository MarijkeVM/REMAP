#### INSTRUCTIONS FOR REMAP (Linux environment) ####

Extract REMAP.rar.
The folder contains a Preprocessing folder, an environment.yml file, examples and the REMAP functions.
 - The environment.yml can be used to create the REMAP conda environment with the REMAP dependencies installed which are used in the REMAP functions (See Installation). 
 - The Preprocessing folder contains functions to create the required input files for the REMAP tool (See Data). 
 - The remaining funcions are the REMAP functions (See Running REMAP and Extract information).

1. Installation:
		Create a conda evironment in which the REMAP functions were developed. The functions are written in Python.2.7.15
		Open a terminal and follow the steps below.
		
		1.a) Download conda: 
             Run the command: (conda will be installed in the current working directory)
				$ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh	
			 
		1.b) Install conda:
			 Run the command:
				$ bash miniconda.sh
		
        1.c) Check if gcc is installed:
			 Run the command:
				$ sudo apt-get install gcc
				
		1.d) Create the evironment:
             Go to the folder containing the REMAP functions
			 Run the command:
			    $ cd Location/of/REMAP_Folder
				$ conda env create -f environment.yml
				
		Activate the REMAP evironment by:
			$ source activate REMAP
			
		Deactivate the REMAP evironment by:
			$ source deactivate
		 

2. Data:

   a) The REMAP function requires a .txt file with a single item (gene) per line with 3 columns separated by tabs:
		column 1 = Manufacturer's ID    
		column 2 = Annotated GeneID's (separated by " /// ")  
		column 3 = Annotated manufacturer's probeset and junction ID's (separated by "|")

      For the HTA-2.0 example, an example is:

		 TC17001298	 SGK494  ///  SPAG5  ///  RP11-192H23.4 	JUC17009694|...|JUC17009742|PSR17017142|...|PSR17017213
		 
      If library files of the corresponding microarray are provided, the required .txt file can be obtained with the provided TCID_GeneID_Probesets.py as shown below.
      The script Merges all information together in the required format illustrated above and creates the required .txt file.
      The necessary files are a transcript annotation file and a probe set annotation file. Further, a prefix for the name of the output file is required.
		
	  USE: python TCID_GeneID_Probesets.py "Location/of/transcript_annotation.csv" "Location/of/probeset_annotation" "location" "prefix"
	   
	  EXAMPLE: $ python TCID_GeneID_Probesets.py "HTA-2_0.na35.2.hg19.transcript.csv" "HTA-2_0.na35.hg19.probeset.csv" "HTA-2_0_Output" "HTA-2_0"
     
      The script is written to work on the structure of the Affymetrix files of the HTA-2.0 microarray but can be easily adapted to different structures. 
     
     
   b) Additional required information is the sequence of the designed probes of the array. We achieve this by combining the manufacturer's .clf and .pgf files with the Probeset_Sequences.py script.
      This creates an out1file which contains the probe sets and their sequences as well as a out2file which contains the line numbers corresponding to a certain probe set for each reading (an indexing file).
	  Further, the number of rows to skip in the .clf file is to be specified. These are the rows starting with "#" and is set to 11 by default.
	  A .pbs file is available for running the script on a supercomputer infrastructure.	
   
      USE: python Probeset_Sequences.py Location/of/clffile Location/of/pgffile skipnrows "location" "out1file" "out2file"
   
      EXAMPLE: $ python Probesets_Sequences.py "HTA-2_0.r3.clf" "HTA-2_0.r3.pgf" 11 "HTA-2_0_Probeset_Sequences.txt" "HTA-2_0_Probeset_SequenceIndices.txt"
           
       
   c) This step is optional: creation of a local copy of the Ensembl database.
      If the database is created locally, go to the Functions/Functions.py and comment line 12 and uncomment line 13.
      A local copy might decrease the computation time.
       
      
    When all information is gathered, REMAP can be run.
	
	
        
3. Running REMAP  (activate the installed REMAP conda environment before running the REMAP function)

	The required input for the REMAP tool is: 
		- The file created by TCID_GeneID_Probesets.py,
		- The out1file created by Probeset_Sequences.py,
		- The out2file created by Probeset_Sequences.py,
		- The name of the output folder in which the produced files are to be saved,
		- The prefix of the output files.
		
	Note the reverse order of out1file and out2file in the command line specifications.	

   a) REMAP.py
   
      This implies running the items one by one on a laptop (not recommended but short examples are included).
      Create an "Output" folder for the output to be saved.
      
      USE: python REMAP_SingleItems.py "TCID_GeneID_Probesets.txt" "Probeset_SequenceIndices.txt" "Probeset_Sequence.txt" "location" "outputfile" "prefix"
      
      EXAMPLE: $ source activate REMAP
	           $ python REMAP.py "HTA-2_0_Examples.txt" "HTA-2_0_Examples_Probeset_SequenceIndices.txt" "HTA-2_0_Examples_Probeset_Sequences.txt" "HTA-2_0_Output" "HTA-2_0"
			   $ python REMAP.py "HJAY_Examples.txt" "HJAY_Examples_Probeset_SequenceIndices.txt" "HJAY_Examples_Probeset_Sequences.txt" "HJAY_Output" "HJAY"
	  
	  The output files consist of: 
	  
	  - "IDS.txt" indicating the TC ID which have run.
	  
	  - "IDS_Completed.txt" indicating the TC IDs which have run and completed.
	  
	  - "GeneInfo.csv" containing the isoform information of the TC IDs as well as the isoform compositions.
	  
	  - "PSRAnnotated.txt" containing the annotation between the PSR IDs and the found annotations to TC IDs and exon IDs for PSR IDs of which all probes match the same exon including the order of the probes and the constructed sequence of the PSR. Multiple annotations are separated by a "|".
	  
	  - "PSRNotAnnotated.txt" containing the PSR IDs without a matching exon ID.
	  
	  - "JUCAnnotated.txt" containing the annotation between the JUC IDs and the found annotations to TC IDs and exon IDs for JUC IDs of which all probes match the same exons. The 5' end exon and 3'end exon are separted by a "|" and multiple annotation are separated by a "//".

	  - "JUCNotAnnotated.txt" containing the JUC IDs without a matching exon IDs.
	  
	  - "GeneAllInfo.txt" containing the TC ID, the PSR ID, the probe ID of the PSR, the exon annotation, the ENSG ID and the provided gene name.
	  
	  - "NotAnnotatedGeneEnsembl.txt" containing the TC IDs of which the provided manufacturer's gene name was unknown to the Ensembl database.
	  
	  - "NotAnnotatedTC.txt" containing the TC IDs without an annotated gene name.
	  

   b) REMAP_remoteDatabase.pbs and REMAP_localDatabase.pbs
      This implies running the items one by one on a supercomputer infrastructure.
      Create an "Output" folder for the output to be saved.
      
	  This requires the installation of the conda environment in the supercomputer infrastructure.
	  
	  Identical output files as for REMAP.py are created.
	  
	  The available example .pbs files have suffixes "_remoteDatabase.pbs" and "_localDatbase.pbs" which are adapated for acessing the Ensembl database remotely and locally respectively.
	  
	
4. Extract information

	Information of the re-annotation can be extracted in a complete format with the script REMAP_ExtractInfo.py. 
	The required input is the GeneAllInfo.txt and GeneInfo.csv files produced by REMAP.py or REMAP_qsub.py. Further, the name of an output file and the prefix of the internally created outputfiles need to be specified.
    A .pbs file is available for running the script on a supercomputer infrastructure.	
	
	 USE: python REMAP_ExtractInfo.py Location/of/GeneAllInfo.txt Location/of/GeneInfo.csv "location" "outputfile" "prefix"
	 
	 EXAMPLE: $ python REMAP_ExtractInfo.py "HTA-2_0_GeneAllInfo.txt" "HTA-2_0_GeneInfo.csv" "HJAY_Output" "HTA-2_0_REMAP.txt" "HTA-2_0"
	 
	 
	 

In case of errors, please feel free to contact marijke.vanmoerbeke@uhasselt.be


      








      
        
   