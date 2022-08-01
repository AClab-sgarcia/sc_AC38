#######################################################################################################################
#Strategies tried to run Seurat V4 in Indar and Nextera. Summary:
#1- normal package installing from R: failed dependencies (see below) in Indar
#2- docker works in Nextera, but in Indar I need root password to initialize the docker.
#3- renv (see other script 0_set_env_indar_scPipeline.R). Need to run this per project. In Indar, issues with openSSL wrote to Asier. Pending
#######################################################################################################################

#------------------------------------------------------------
#Trial1: normal package installing from R in Indar:
#------------------------------------------------------------
##################################################
#ask for iterative job in Indar
##################################################
salloc -N 1 -n 1 --mem=30G -t 4:00:00 --partition=FAST --job-name=interactive
ssh cn04
squeue

##################################################
#execute R v4 in indar from the folder where the renv environment is located
#initalize env 
##################################################
/share/apps/R/R-4.1.2/bin/R

##################################################
#install packages and dependencies. Could not make it
##################################################
install.packages("Seurat",dependecies=T)
ERROR: dependencies ‘cowplot’, ‘future’, ‘future.apply’, ‘ggrepel’, ‘ggridges’, ‘httr’, ‘igraph’, ‘leiden’, ‘miniUI’, ‘patchwork’, ‘plotly’, ‘png’, ‘RcppAnnoy’, ‘reticulate’, ‘ROCR’, ‘Rtsne’, ‘scattermore’, ‘sctransform’, ‘SeuratObject’, ‘uwot’, ‘RcppEigen’ are not available for package ‘Seurat’

suppressPackageStartupMessages(require("dplyr"))
suppressPackageStartupMessages(require("scater"))
suppressPackageStartupMessages(require("Seurat"))
suppressPackageStartupMessages(require("patchwork"))
suppressPackageStartupMessages(require("cowplot"))
suppressPackageStartupMessages(require("ggplot2"))
suppressPackageStartupMessages(require("celldex"))
suppressPackageStartupMessages(require("purrr"))
suppressPackageStartupMessages(require("SingleR"))
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
suppressPackageStartupMessages(require("DoubletFinder"))
suppressPackageStartupMessages(require("scales"))
suppressPackageStartupMessages(require("biomaRt"))
suppressPackageStartupMessages(require("SeuratWrappers"))
# remotes::install_github("satijalab/seurat-wrappers")
suppressPackageStartupMessages(require("biomaRt"))
suppressPackageStartupMessages(require("SeuratDisk"))
# remotes::install_github("mojaveazure/seurat-disk")
source("utils/custom_seurat_functions.R")

#------------------------------------------------------------
#Trial 2: Docker in Nextera
#------------------------------------------------------------
##################################################
# docker https://satijalab.org/seurat/articles/install.html
##################################################
#systemctl start docker
#docker pull satijalab/seurat:latest
#To use as a base image in a new Dockerfile:
#FROM satijalab/seurat:latest
docker images
docker run -ti -v /vols/GPArkaitz_bigdata/sgarcia/Single-cell_AC38_AC50/:/docker/ --rm docker.io/satijalab/seurat R
# we need to setwd "docker" as we have mounted the docker there with "/:/docker/" when running the docker in the previous code line
setwd("docker")
source("./0_Scripts/2_scPipeline_Integration_13052022.R")

#------------------------------------------------------------
#Trial 3: renv in Indar and Nextera
#------------------------------------------------------------
##################################################
# 3.1 ask for iterative job in Indar
##################################################
salloc -N 1 -n 1 --mem=30G -t 4:00:00 --partition=FAST --job-name=interactive
#cn04 
ssh cn04
squeue

##################################################
#execute R v4 in indar from the folder where the renv environment is located
#initalize env 
##################################################
cd /vols/GPArkaitz_bigdata/sgarcia/Single-cell_AC38_AC50_/
/share/apps/R/R-4.1.2/bin/R
renv::init()
#1
q()

#Now install the packages again in linux (only nee to be done ones)
/share/apps/R/R-4.1.2/bin/R
renv::install("chris-mcginnis-ucsf/DoubletFinder")
renv::install("scales")
renv::install("cowplot")
renv::install("purrr")
renv::install("ggplot2")
renv::install("scCATCH")
renv::install("patchwork")
renv::install("data.table")
renv::install("dplyr")
install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("celldex")
BiocManager::install("MAST")
BiocManager::install("scater")

renv::install("hadley/devtools")
renv::install("r-lib/textshaping")

renv::install("Seurat")
#|||||||||||||||||||||||||||||||||||||||||||||||||||
Issues with openSSL wrote to Asier. Pending
#|||||||||||||||||||||||||||||||||||||||||||||||||||
#did not try these below:
renv::install("mojaveazure/seurat-disk")
renv::install("satijalab/seurat-wrappers")

#last step:
renv::snapshot() 

##################################################
# 3.2 renv in Nextera. Issues with RccpTOML
##################################################
cd /vols/GPArkaitz_bigdata/sgarcia/Single-cell_AC38_AC50/
/opt/R/R-4.1.2/bin/R
renv::init()
#1
q()
