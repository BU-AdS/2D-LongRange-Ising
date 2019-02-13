#!/bin/bash -l

#$ -l h_rt=120:00:00
#$ -l mem_total=125G
#$ -m a
#$ -j y
#$ -N J=0.4406867935

module load gsl
J=0.4406867935/template_ising.sh

