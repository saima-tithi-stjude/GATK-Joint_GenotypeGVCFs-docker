#!/usr/bin/env nextflow

/***********************************************************************************************
 This Nextflow Pipeline as the second part, which can be called to generate the cohort joint VCF
 Input : a list for gVCFs per sample
 Output: cohort VCF 
 Written by: Wenchao Zhang
 Center for Applied Bioinformatics, St. Jude Children Research Hosptial
 Date: 01/11/2022
 Updated by Saima Sultana Tithi, added option for running it in a docker/singularity container, 11/18/2024
***********************************************************************************************/


REFERENCE=""
REFERENCE_PATH=""
REFERENCE_BASE=""
KNOWN_SITE1=""
KNOWN_SITE2="" 
KNOWN_SITE3="" 

EXOME_BED="" 
INTERVALS=""
withChr=true
PICARD="/gatk/picard.jar"

hapmap=""
omni=""
thousandgenome_snp=""
dbSNP=""
mills_indel=""
axiom=""
ICAnnotation=""
DPAnnotation=""
MQAnnotation="-an MQ"

def parse_config_parameters() {
// Parse Nextflow Configurations and Paramters 
    if( params.genome_build == 'hg38' )
    {  
        REFERENCE = params.resourceFiles + "/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
        KNOWN_SITE1 = params.resourceFiles + "/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
        KNOWN_SITE2 = params.resourceFiles + "/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KNOWN_SITE3 = params.resourceFiles + "/hg38/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
        EXOME_BED = params.resourceFiles + "/hg38/sureSelectV6_hg38_regions.bed"
        INTERVALS = params.resourceFiles + "/hg38/hg38_intervals.txt"

        hapmap = params.resourceFiles + "/hg38/GRCh38_no_alt/GATK/hapmap_3.3.hg38.vcf.gz"
        omni   = params.resourceFiles + "/hg38/GRCh38_no_alt/GATK/1000G_omni2.5.hg38.vcf.gz"
        thousandgenome_snp = params.resourceFiles + "/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        dbSNP = params.resourceFiles + "/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"

        mills_indel = params.resourceFiles + "/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        axiom = params.resourceFiles + "/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
    }
    else if( params.genome_build == 'hg19' )
    {
        REFERENCE = params.Build_hg19.REFERENCE
        KNOWN_SITE1 = params.Build_hg19.KNOWN_SITE1
        KNOWN_SITE2 = params.Build_hg19.KNOWN_SITE2
        KNOWN_SITE3 = params.Build_hg19.KNOWN_SITE3
        EXOME_BED = params.Build_hg19.EXOME_BED
        INTERVALS = params.Build_hg19.INTERVALS

        hapmap = params.Build_hg19.hapmap 
        omni   = params.Build_hg19.omni
        thousandgenome_snp = params.Build_hg19.thousandgenome_snp
        dbSNP = params.Build_hg19.dbSNP
        
        mills_indel=params.Build_hg19.mills_indel
        axiom=params.Build_hg19.axiom   
    }
    else if( params.genome_build == 'b37' ) 
    {
        withChr=false
        REFERENCE = params.Build_b37.REFERENCE
        KNOWN_SITE1 = params.Build_b37.KNOWN_SITE1
        KNOWN_SITE2 = params.Build_b37.KNOWN_SITE2
        KNOWN_SITE3 = params.Build_b37.KNOWN_SITE3
        EXOME_BED = params.Build_b37.EXOME_BED
        INTERVALS = params.Build_b37.INTERVALS
        
        hapmap = params.Build_b37.hapmap 
        omni   = params.Build_b37.omni
        thousandgenome_snp = params.Build_b37.thousandgenome_snp
        dbSNP = params.Build_b37.dbSNP

        mills_indel=params.Build_b37.mills_indel
        axiom=params.Build_b37.axiom   
    }
    else
    {
        error "Invalid geome version: ${params.genome_build}"
        exit 1
    }

    if( params.popType == "POP" && params.ASOption != "AS")
    { 
        ICAnnotation="-an InbreedingCoeff"
    }

    if ( params.data_type == "WHOLE_GENOME" && params.ASOption != "AS" )
    {
        DPAnnotation="-an DP"
    }

    REFERENCE_PATH = file(REFERENCE).getParent()
    REFERENCE_BASE = file(REFERENCE).getBaseName()
}

def DispConfig() {
 log.info """
Welocme to run Nextflow Pipeline GenotypeGVCFs.nf  
Your configuration are the following:
  GVCF_filelist    : ${params.sampleList}
  genome_build     : ${params.genome_build}
  data_type        : ${params.data_type}
  vqsr_recal       : ${params.vqsr_recal}
  ASOption         : ${params.ASOption}
  popType          : ${params.popType}
  VQSR_indelTranche: ${params.VQSR_indelTranche}
  VQSR_snpTranche  : ${params.VQSR_snpTranche}

The parsing information based your configuration are the following:
  KNOWN_SITE1      : $KNOWN_SITE1
  KNOWN_SITE2      : $KNOWN_SITE2     
  KNOWN_SITE3      : $KNOWN_SITE3
  EXOME_BED        : $EXOME_BED
  INTERVALS        : $INTERVALS   
  hapmap           : $hapmap 
  omni             : $omni
  thousandgenome_snp : $thousandgenome_snp
  dbSNP            : $dbSNP
  mills_indel      : $mills_indel
  axiom            : $axiom 
  withChr          : $withChr
 """
 // exit 0    
}

def helpMessage() {
  log.info """
        Welocme to run Nextflow Pipeline GenotypeGVCFs.nf 
        Usage:
        A typical command for running the pipeline is as follows:
        nextflow run GenotypeGVCFs.nf -profile cluster -w ~/Nextflow_work --outdir ./
        A command with more configurable arguments can be           
        nextflow run GenotypeGVCFs.nf -profile cluster -w ~/Nextflow_work --GVCF_filelist ./nova_GVCF_list.txt --project Ptest --genome_build hg38 --data_type WHOLE_GENOME vqsr_recal Y --outdir ./

        Configurable arguments:
        --sampleList                     GVCF file list that you want to do the joint cohort calling analysis                
        --genome_build                   Reference genome version. hg38(default) | hg19 | b37 
        --data_type                      WHOLE_GENOME(default) | EXOME
        --vqsr_recal                     VQSR or not. Y(default) | N
        --ASOption                       Allele Specific Option. AS(Default) | (Empty)
        --popType                        Population datatype, InbreedingCoeff used or not. notPOP(Default) | POP
        --VQSR_indelTranche              VQSR INDEL Tranche level for truth sensitivity filter
        --VQSR_snpTranche                VQSR SNP Tranche level for truth sensitivity filter
        --help | h                       This usage statement.
       
        Note:
         All arguments can be configured in the command line interface or specified by the nextflow.config (default)  
        """
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

//Parse the input Parameters and configs. 
parse_config_parameters() 
// Display the configuration
DispConfig() 

//
// Convert the input Parameters/argument into the correct format 
//
gvcf_in_ch = Channel            
			   .fromPath(params.sampleList)               
               .splitText()
               .map{it.replaceFirst(/\n/, "")}
               .map{ file(params.gvcfFolder + "/" + it) }   //map the file path string into file object, then can extract the file information.  

INTERVAL_in_ch = Channel
                 .fromPath(INTERVALS)
                 .splitText()
                 .map{it.replaceFirst(/\n/, "")}
                 .map{it -> tuple( it.split(":")[0], it, it.replaceFirst(/:/, "-"))}

/***************************************************************************************************************************************
* This process is responsible for generating genome reference indexing
***************************************************************************************************************************************/
process genome_indexing {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    
    input:
    val genome from REFERENCE

    output:
    set file("Ref.fa"), file("Ref.fa.fai"), file("Ref.dict") into (Ref_faidx_corhort_ch, Ref_faidx_vqsr_snp_ch, Ref_faidx_vqsr_indel_ch, Ref_faidx_vqsr_apply_ch)
    
    script:
    if( file("${genome}.fai").isEmpty() || file("${REFERENCE_PATH}/${REFERENCE_BASE}.dict").isEmpty())
    """
      ln -s ${genome} Ref.fa
      samtools faidx Ref.fa
      java -jar ${PICARD} \
      CreateSequenceDictionary \
      R= "Ref.fa" \
      O= "Ref.dict" 
    """
    else  
    """
      ln -s ${genome} Ref.fa
      ln -s ${genome}.fai Ref.fa.fai
      ln -s ${REFERENCE_PATH}/${REFERENCE_BASE}.dict Ref.dict
    """
}

/**************************************************************************************************************************************
* This process is responsible for generating the gVCF indexing file (.tbi) if missing   
*************************************************************************************************************************************/
process gvcf_indexing {    
    container 'stithi/gatk_jointgenotype_gvcfs:v1'  

    input:      
    set filename, parent, gvcf_file from gvcf_in_ch.map{ file -> tuple(file.name, file.parent, file) }  
    
    output:  
    file("${filename}") into gvcf_ch   
    file("${filename}.tbi") into gvcf_tbi_ch 
    
    script:
    if (file("${parent}/${filename}.tbi").isEmpty())
    """
    ln -s ${parent}/${filename} ${filename} 
    bcftools index ${parent}/${filename} -t
    """
    else
    """
    ln -s ${parent}/${filename} ${filename}
    ln -s ${parent}/${filename}.tbi ${filename}.tbi
    """
}

/***************************************************************************************************************************************
* This process is responsible for collecting all gVCFs and generating the 
interval level Genomic Database , such as test-project.chr21-42038964-43206713.db
***************************************************************************************************************************************/
process GenomicsDBImport {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    //publishDir "${params.outdir}/${params.project}/GenomicsDBImport", mode: 'symlink', overwrite: true
    // time { (24.hour + (24.hour * task.attempt)) } // First attempt 24h, second 48h, etc
    memory { (16.GB * task.attempt) } // First attempt 32GB, second 64GB, etc
    errorStrategy 'retry'
    maxRetries 3

    input:	
    file(gvcf) from gvcf_ch.collect()
    file(gvcf_tbi) from gvcf_tbi_ch.collect()
    set val(Chr), val(INTERVAL), val(INTNAME) from INTERVAL_in_ch
   
	output:
    set Chr, INTERVAL, INTNAME, file ("${params.project}.${INTNAME}.db") into gendb_ch		
  
    script:
    """
    gatk --java-options "-Xmx24g -Xms24g" \
    GenomicsDBImport \
    ${gvcf.collect { "-V $it " }.join()} \
    -genomicsdb-workspace-path ${params.project}.${INTNAME}.db \
    -L ${INTERVAL} \
    --reader-threads 5 \
    --batch-size 50
    """
}	

/*********************************************************************************************************************************
This process is responsible for launching GenotypeGVCFs on the previously created interval level genomic Database and output the 
interval level joint called VCF. 
**********************************************************************************************************************************/
process GenotypeGVCFs {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
	//publishDir "${params.outdir}/${params.project}/GenotypeGVCF", mode: 'symlink', overwrite: true
    // time { (48.hour + (48.hour * task.attempt)) } // First attempt 48h, second 96h, etc
    memory { (16.GB * task.attempt) } // First attempt 32GB, second 64GB, etc
    errorStrategy 'retry'
    maxRetries 3

    input:
	set val(Chr), val(INTERVAL), val(INTNAME), file(INTERVAL_Genomic_DB) from gendb_ch
   	set file(genome), file(faidx), file(fadict) from Ref_faidx_corhort_ch
    val DBsnp from KNOWN_SITE2

	output:
    set Chr, file("${params.project}.${INTNAME}.cohort.vcf.gz"), file("${params.project}.${INTNAME}.cohort.vcf.gz.tbi") into cohort_interval_vcf_ch

    script:
	"""
    gatk --java-options "-Xmx24g -Xms24g" \
    GenotypeGVCFs \
    -R ${genome} \
    -V gendb://${INTERVAL_Genomic_DB} \
    -O ${params.project}.${INTNAME}.cohort.vcf.gz \
    -L ${INTERVAL} \
    -G StandardAnnotation \
    -G AS_StandardAnnotation \
    -D $DBsnp \
    --only-output-calls-starting-in-intervals  
	"""
}	

/******************************************************************************************************************************************
/* This Process is responsible for Merging the interval VCFs (sub chrosome level, e.g. test-project.chr21-42038964-43206713.cohort.vcf.gz) 
 into chrosome level VCF. 
********************************************************************************************************************************************/
process Merge_Intervals_ToChromVCF {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.project}/Merge_Intervals_ToChromVCF", mode: 'copy'
    input:
	set val(Chr), file (vcfs), file (vcftbis) from cohort_interval_vcf_ch.groupTuple()
       
	output:
    set Chr, file("${params.project}.${Chr}.vcf.gz"), file("${params.project}.${Chr}.vcf.gz.tbi") into (cohort_chr_vcf_ch, vqsr_chr_vcf_ch, chr_vcf_recal_ch)

    script:    
    """
    gatk --java-options '-Xmx24g -Xms24g' \
    MergeVcfs \
    ${vcfs.collect{"--INPUT $it " }.join()} \
    -O ${params.project}.${Chr}.vcf.gz
    """
}
