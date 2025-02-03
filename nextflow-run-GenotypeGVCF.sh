module load nextflow/21.10.5
nextflow run GenotypeGVCFs-docker.nf -profile cluster_singularity -w /scratch_space/stithi/Nextflow_work/GMKF-joint-calling -c nextflow.config -with-report nextflow-report.html -with-timeline nextflow-timeline.html -with-dag nextflow-dag.png
