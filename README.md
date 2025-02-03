# A nextflow based pipeline to generate the cohort joint VCF using Singularity container #
module load nextflow/21.10.5
nextflow run GenotypeGVCFs-docker.nf -profile cluster_singularity -w /scratch_space/stithi/Nextflow_work/testing/ -c nextflow.config -with-report nextflow-report.html -with-timeline nextflow-timeline.html -with-dag nextflow-dag.png

Notes to nextflow command line configuration options
1.) Nextflow related configuration options (https://www.nextflow.io/docs/latest/config.html ):
-profile cluster_singularity  :     select cluster_singularity for running in a singularity container
-w  /scratch_space/stithi/Nextflow_work/ :   the directory for caching nextflow's process result.    
-resume: nextflow will resume the pipeline after you fix the error/bugs.
-c nextflow.config :  use my specific configuration file nextflow.config   

2.) Pipeline related parameter options:
--h:  Help Message
--help:  Help Message
--GVCF_filelist:   gVCF filelist for GenotypeGVCFs.nf's input
--outdir:    Pipeline result output directory 
--project:  Correspond to one pipeline's run.        
--genome_build:  hg38|hg19|b37
--data_type:    WHOLE_GENOME|EXOME
--vqsr_recal:  Y/N 
--VQSR_indelTranche:    VQSR INDEL Tranche  
--VQSR_snpTranche:   VQSR SNP Tranche  

