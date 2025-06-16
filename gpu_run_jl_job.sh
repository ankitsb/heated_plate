#!/usr/bin/env bash
#PBS -N heat_plat
#PBS -P ds2
#PBS -q gpuvolta
#PBS -l walltime=48:00:00
#PBS -l ncpus=12
#PBS -l ngpus=1
#PBS -l mem=30GB
#PBS -l jobfs=100GB
#PBS -l storage=scratch/ds2+gdata/hh5+gdata/jk72
#PBS -j oe
#PBS -l wd

module purge
module load julia/1.9.1
module load cuda

julia heated_plat_code.jl >out.dat 2>&1
