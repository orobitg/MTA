#SGE arguments
#$ -V
#$ -N Name
#$ -S /bin/bash
#$ -o /home/user/out.txt
#$ -e /home/user/error.txt
#$ -pe mpi 16 #Environtment and number of cpus
#$ -M mail@mail.cat
#$ -m e
#$ -v MPIR_HOME=/home/user/mpich3


mpirun -n 16 bin/mta.mpi bin/mta.mpi -seq /home/user/seq.fasta -ntree 100 -msa t_coffee -score sp -output /home/user/results -mpi
