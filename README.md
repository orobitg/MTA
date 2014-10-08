MTA-NF
======

MTA is a Multiple Sequence Alignment pipeline to align multiple guide-tree variations for the same input sequences, evaluating the alignments obtained and selecting the best one as result. MTA method can be applied to any progressive method that accepts guide trees as an input parameter (in this version accepts T-Coffee, ClustalW and Clustal-Omega). In addition, the method also allows different evaluation metrics to select the best multiple guide trees (sp, normd). The aim is to find a variation of the original tree that provides a more accurate alignment than the original one produced. MTA is implemented in C, and there is a parallel version usign MPI libraries.


Quick start 
-----------

Clone the git repository on your computer with the following command:

    $ git clone git@github.com:orobitg/MTA.git
    
Run the script install.sh to download, compile and install the tools listed below (dependecies). The install.sh script also compiles MTA and creates two binaries in the bin folder:
	
	+ mta: The serial version of MTA.
	+ mta.mpi: The mpi version of MTA.

When done, move in the project root folder named `MTA`, 
which contains an example dataset in the `tutorial` folder. 

Launch the pipeline by entering the following command 
on your shell terminal:

    $ bin/mta -seq ./tutorial/sample.fa

To re-compile the mta and mta.mpi binaries:
	
	`make -f Makefile`
	`make -f Makefile.mpi`
    
By default the pipeline is executed against the provided tutorial dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program command line.

Pipeline parameters
-------------------

**-seq**  
   
* The location of the sequences fasta file.
  	`$ bin/mta -seq /home/user/seq/example.fa`

**-ntree**  
   
* Number of of guide trees to align and evaluate.
	* By default is set to 10 guide trees  
  	`$ bin/mta -seq /home/user/seq/example.fa -ntree 100`

**-msa**  
   
* The MSA tool to produce the alignemnts.
* Two options: t_coffee, clustalw, clustalo.
	* By default is set to 'T-Coffee'  
	`$ bin/mta -seq /home/user/seq/example.fa -msa t_coffee`

**-score**  
   
* The evaluation score to choose an alignemnt.
* Two options: sp, normd, tcs (t_coffee score)
	* By default is set to 'Sum-of-Pairs (sp)'    
  	`	$ bin/mta -seq /home/user/seq/example.fa -score sp`

**-gop** 
   
* Sets the Gap Opening Penalty. Only used with SP score.  
	* By default is set to -11    
  	`$ bin/mta -seq /home/user/seq/example.fa -score sp -gop -11`

**-gep** 
   
* Sets the Gap Extended Penalty. Only used with SP score.  
	* By default is set to -1   
  	`$ bin/mta -seq /home/user/seq/example.fa -score sp -gep -1`

**-matrix** 
   
* Sets the Distance Matrix. Only used with SP score.  
* Options: blosum30mt, blosum40mt, blosum45mt, blosum50mt, blosum55mt, blosum62mt, blosum80mt, idmat, dna_idmat, pam120mt, pam160mt, pam250mt, pam350mt, md_40mt, md_120mt, md_250mt, md_350mt.
	* By default is set to 'blosum62mt'    
  	`$ .bin/mta -seq /home/user/seq/example.fa -score sp -matrix blosum62mt`  
  
**-output** 
   
* Specifies the folder where the results will be stored for the user.
* It does not matter if the folder does not exist.
  	* By default is set to MTA-NF's folder: './'  
  	`$ bin/mta -seq /home/user/seq/example.fa -output /home/user/my_results`

**-mpi**

*You must specify this argument if you are running the parallel version.
	`$ mpirun -np 4 bin/mta.mpi -seq /home/user/seq/example.fa -mpi`

Run locally 
------------

**Command line**

Install all the dependencies running the bash script install.sh.

	`$ bash install.sh`

You can re-compile the binarie doing:

	`$ make -f Makefile`

Run MTA command line indication all the input parameters. For example:

	$ bin/mta -seq /home/user/seq.fasta -ntree 100 -msa t_coffee -score sp -output /home/user/results

  
Cluster support - MPI version
-----------------------------

MTA also supports MPI libraries. Thus it is possible to execute it on your computer or any cluster with mpicc compiler and mpirun command installed.

Install all the dependencies running the bash script install.sh.

	$ bash install.sh`

You can re-compile the binarie doing:

	`$ make -f Makefile.mpi`

Run the mpi command line indication all the input parameters (run locally with mpirun). For example:

	$ mpirun -np 4 bin/mta.mpi -seq /home/user/seq.fasta -ntree 100 -msa t_coffee -score sp -output /home/user/results -mpi

To submit in a SGE/PBS cluster create a script similar to the script found in tutorial/sge_script.sh. Run the script through the `qsub` command.


Dependencies 
------------

 * T-Coffee - http://www.tcoffee.org/Projects/tcoffee/index.html
 * ClustalW - http://www.clustal.org/clustal2/
 * Normd - ftp://ftp-igbmc.u-strasbg.fr/pub/NORMD/
 * Clustalo - http://www.clustal.org/omega/
 * MPICH - http://www.mpich.org/downloads	

