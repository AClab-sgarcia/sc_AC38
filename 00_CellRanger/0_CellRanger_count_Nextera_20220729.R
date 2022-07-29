##################################################
#cellranger VERSION INFO
##################################################
#cell ranger version in Nextera is latest in July 2022: cellranger v7.0.0 . In this version the estimated cell number is not longer needed and it is estimated automatically.
#cell ranger version in Indar is V6.1.0
# refernce genomes are from 2020 (latest) according to https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build
##################################################

##################################################
#including intron reads with option?
##################################################
# From https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjL27jJsYb2AhWRZMAKHRSIDPgQFnoECAQQAQ&url=https%3A%2F%2Fassets.ctfassets.net%2Fan68im79xiti%2FawNZTarmwqmwxKcvDv9wv%2Fb05500661e36290b8ce59689ff889ea8%2FCG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf&usg=AOvVaw11gnGq4fzMOsTX_h09r4nU
#Counting these intronic reads alongside exonic reads can yield significant increases in assay sensitivity, especially for nuclei samples. Moreover, intronic UMI inclusion leads to the detection of genes that have no exonic UMIs in any cell. On the other hand, this internal priming appears to result in the assignment of multiple UMIs to single transcripts. This causes a gene length bias for intronic UMIs that should be taken into account when considering including intronic UMIs in Single Cell Gene Expression analysis.

##################################################
# HOW to run in Nextera:
# 1- input folder in DATA_SHARED
# 2- output folder inside each users folder, the script will generate the outputs within "CellRangerCount" folder. 
# WARNING!! IF the output file is same as input file, the RAW data will be overwritten. 
# 3- This script will generate a file listing all commands (one per sample) to send it to the home-made job scheduler for Nextera (homeMadeParallel.sh)
##################################################

/opt/R/R-4.1.2/bin/R

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
project_name <- "AC38"
samples <- c("S_01", "S_02")
dir_outfiles <- "/vols/GPArkaitz_bigdata/sgarcia/sc_AC38/" # must be user's path
dir_infiles <- "/vols/GPArkaitz_bigdata/DATA_shared/AC-38_10x-3RNA/FASTQs/" # must be DATA_shared path

# Organism (options: "human" || "mouse")
organism <- "mouse"

if (organism == "human") {
  genome <- "refdata-gex-GRCh38-2020-A"
} else if (organism == "mouse"){
  genome <- "refdata-gex-mm10-2020-A"
}
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
dir.create(file.path(dir_outfiles,"CellRangerCount"))
dir_outfiles <- paste(dir_outfiles,"CellRangerCount",sep='')
setwd(dir_outfiles)

system("cp /vols/GPArkaitz_bigdata/DATA_shared/utils/homeMadeParallel.sh ./")
dir_genome <- "/vols/GPArkaitz_bigdata/DATA_shared/Genomes/"

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
memory <- c("48")
cpu <- 8
numb_process <- 3 # how many processes will be run in parallel. Take into account how much memory each job takes (total cpus in Nextera are 32 and 188GB of memory)
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

filename <- paste(project_name,".sh",sep='')
#since the loop appends to the previously existing file, make sure to remove any previous file with the same name as done by this  line
system(paste("rm ",filename,sep=' '))
for (f in 1:length(samples)) {
  command <- paste("cellranger count --id=",samples[f]," --transcriptome=",dir_genome,genome," --fastqs=",dir_infiles,"/"," --sample=",samples[f]," --localcores=",cpu," --localmem=",memory," --no-bam",sep='') 
  cat(
    c(command),
    file=filename,sep = "\n",append=T);
  #system(paste("chmod 777",filename,sep=' '))
}
system(paste("./homeMadeParallel.sh",filename,numb_process,sep=' '))

#####################################################################
#cd /storage/AClab/10_AMD1_splicing_CRG/vast-tools/
#./homeMadeParallel.sh 1_align_AC27_AC19_CRG.sh 3
#####################################################################
