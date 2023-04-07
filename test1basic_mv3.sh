#!/bin/bash
#SBATCH --job-name=LMPPHN         # Job name
#SBATCH -N 5 #Number of nodes
#SBATCH --ntasks-per-node=48  #Number of core per node
#SBATCH --error=job.%J.err  #Name of output file
#SBATCH --output=job.%J.out #Name of error file

#SBATCH --partition=medium #specifies queue name(standard is the default partition if you does not specify any partition job will be submitted using default partition) other partitions You can specify hm and gpu


module load spack
source /home/apps/spack/share/spack/setup-env.sh
spack load lammps@20220623/d6ambcm
module load gcc-8.3.0-gcc-11.2.0-lyhn5yb


#This is for normal lammps
#spack load lammps /ajb3sh4
###################################
#This is for phonon package in lammps
#spack load lammps@20220623/wmkb6ui
#module load gcc-8.3.0-gcc-11.2.0-lyhn5yb
####################
##Did not work##
#ml intel/2020
#ml intel/oneapi/mpi/2021.2.0 
#ml gcc/10.2
###################
#source /opt/ohpc/pub/compiler/intel/2020/compilers_and_libraries/linux/mkl/bin/mklvars.sh intel64

#spack load lammps@20220623/sobl24j
#spack load gcc/bwq7xaa
#spack load openmpi/7cvclvr

#export PATH="/home/a.deshmukh/apps_AD/lammps-3Mar20/src/:$PATH"
#export PATH="/home/param.singh/apps/LAMMPS/lammps-29Oct20/build:$PATH"
export OMP_NUM_THREADS=5

mpirun -np 240 lmp -in in.lammps

