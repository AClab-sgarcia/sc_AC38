# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct

# cellranger count --help

# cellranger-count
# Count gene expression (targeted or whole-transcriptome) and/or feature barcode reads from a single sample and GEM well
 
# USAGE:
#     cellranger count [FLAGS] [OPTIONS] --id  --transcriptome 
 
# FLAGS:
#         --no-bam                  Do not generate a bam file
#         --nosecondary             Disable secondary analysis, e.g. clustering. Optional
#         --include-introns         Include intronic reads in count
#         --no-libraries            Proceed with processing using a --feature-ref but no Feature Barcode libraries
#                                   specified with the 'libraries' flag
#         --no-target-umi-filter    Turn off the target UMI filtering subpipeline. Only applies when --target-panel is
#                                   used
#         --dry                     Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
#         --disable-ui              Do not serve the web UI
#         --noexit                  Keep web UI running after pipestance completes or fails
#         --nopreflight             Skip preflight checks
#     -h, --help                    Prints help information

# Run cellranger count

# To run cellranger count, you need to specify an --id. This can be any string, which is a sequence of alpha-numeric characters, underscores, or dashes and no spaces, that is less than 64 characters. 
# Cell Ranger creates an output directory that is named using this id. This directory is called a "pipeline instance" or pipestance for short.

# The --fastqs should be a path to the directory containing the FASTQ files. If you demultiplexed your data using cellranger mkfastq, you can use the path to fastq_path directory in the outs from the pipeline. 
# If there is more than one sample in the FASTQ directory, use the --sample argument to specify which samples to use. This --sample argument works off of the sample id at the beginning of the FASTQ file name. 
# It is unnecessary for this tutorial run because all of the FASTQ files are from the same sample, but it is included as an example. The last argument needed is the path to the --transcriptome reference package. 
# Be sure to edit the file paths in red in the command below.

# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/7.0/using/count

cellranger count 
    --id
    --fastqs
    --libraries (hen using this argument, --fastqs and --sample must not be used. This argument should not be used when performing gene expression-only analysis; use --fastqs instead.)
    --sample
    --transcriptome
    --feature-ref
    --target-panel
    --no-target-umi-filter (Optional)
    --expect-cells (Optional) (Override the pipeline’s auto-estimate.)
    --force-cells (Optional)
    --include-introns (Optional)
    --nosecondary (Optional)
    --no-bam (Optional) 
    --no-libraries (Optional) 
    --chemistry (Optional) (NOTE: by default the assay configuration is detected automatically, which is the recommended mode. You should only specify chemistry if there is an error in automatic detection)
    --r1-length (Optional) 
    --r2-length (Optional) 
    --lanes (Optional) 
    --localcores (Optional) (Restricts cellranger to use specified number of cores to execute pipeline stages. By default, cellranger will use all of the cores available on your system)
    --localmem (Optional) (Restricts cellranger to use specified amount of memory (in GB) to execute pipeline stages. By default, cellranger will use 90% of the memory available on your system())
    --check-library-compatibility (Optional) 

# When the output of the cellranger count command says, “Pipestance completed successfully!”, this means the job is do

cellranger count 
    --id=
    --fastqs=/vols/GPArkaitz_bigdata/DATA_shared/AC-56_scRNAseq/
    --sample
    --transcriptome=/vols/GPArkaitz_bigdata/DATA_shared/Genomes/
    --localcores=8 
    --localmem=64
    --no-bam 


