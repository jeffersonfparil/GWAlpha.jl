#!/bin/bash

###################################################################
### PREPARE SPARTAN FOR QUANTINEMO SIMULATION AND DATA ANALYSIS ###
###################################################################

### SAMPLE USAGE
# projectID=punim0543
# username=$USER
# DIR=/data/cephfs/${projectID}/${username}
# # DIR=/data/Lolium/Quantitative_Genetics/LOLSIM_2019
# ./GPWASim_00_installation.sh $DIR

### INPUT
DIR=$1

### create installation directory
mkdir ${DIR}
mkdir ${DIR}/Softwares
mkdir ${DIR}/Scripts
mkdir ${DIR}/Output
echo -e "Your installation directory is in $DIR"

### prepare quantiNemo2, genomic_prediction.git, popoolation2, and julia1 libraries
###### download quantiNemo2 if not yet installed
if [ $(ls ${DIR}/Softwares/quantinemo_linux | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading QuantiNemo2"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://www2.unil.ch/popgen/softwares/quantinemo/files/quantinemo_linux.zip
  unzip quantinemo_linux.zip
  rm quantinemo_linux.zip
  echo "########################################"
fi
###### clone the genomic_prediction git repo not yet installed
if [ $(ls ${DIR}/Softwares/genomic_prediction | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Cloning genomic_prediction.git repository"
  echo "########################################"
  cd ${DIR}/Softwares
  git clone https://gitlab.com/jeffersonfparil/genomic_prediction.git
  echo "########################################"
fi
###### download popoolation2 if not yet installed
if [ $(ls ${DIR}/Softwares/popoolation2_1201 | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading Popoolation2"
  echo "########################################"
  cd ${DIR}/Softwares
  wget https://sourceforge.net/projects/popoolation2/files/popoolation2_1201.zip
  unzip popoolation2_1201.zip -d popoolation2_1201
  rm popoolation2_1201.zip
  echo "########################################"
fi
###### download NPSTAT if not yet installed
if [ $(ls ${DIR}/Softwares/npstat | wc -l) -eq 0 ]
then
  echo "########################################"
  echo "Downloading NPSTAT"
  echo "########################################"
  cd ${DIR}/Softwares
  git clone https://github.com/lucaferretti/npstat.git
  cd npstat/
  make || module load GSL/2.5-intel-2018.u4; make
  echo "########################################"
fi
# ###### download julia 1.1 if not yet installed
# if [ $(ls ${DIR}/Softwares/julia-1.1.1 | wc -l) -eq 0 ]
# then
#   echo "########################################"
#   echo "Downloading Julia 1.1.1"
#   echo "########################################"
#   cd ${DIR}/Softwares
#   wget https://julialang-s3.julialang.org/bin/linux/x64/1.1/julia-1.1.1-linux-x86_64.tar.gz
#   tar -xvzf julia-1.1.1-linux-x86_64.tar.gz
#   rm julia-1.1.1-linux-x86_64.tar.gz
#   echo "########################################"
# fi
##### NOTE: INSTALL R PACKAGES MANUALLY!!! (FOR julia's RCall and ggmix)
# module load R/3.5.2-GCC-6.2.0
# R
# install.packages("pacman")
# pacman::p_load_gh('sahirbhatnagar/ggmix')
##### NOTE: INSTALL PYTHON'S MATPLOTLIB MANUALLY!!! (FOR julia's PyCall and Plots' plotting)
# module load Python/2.7.13-GCC-6.2.0-bare
# pip install matplotlib --user
#### NOTE: INSTALL JULIA LIBRARIES MANUALLY
# module load R/3.5.2-GCC-6.2.0
# module load Julia/1.1.1-spartan_gcc-6.2.0.lua
# julia
# using Pkg
# Pkg.add([
# "CSV",
# "Distributions",
# "DataFrames",
# "DelimitedFiles",
# "LinearAlgebra",
# "Optim",
# "Plots",
# "PyPlot",
# "UnicodePlots",
# "ColorBrewer",
# "ProgressMeter",
# "Statistics",
# "StatsBase",
# "MultipleTesting",
# "SortingAlgorithms",
# "PyCall",
# "Conda",
# "GLM",
# "Lasso",
# "GLMNet",
# "RCall"
# ])
