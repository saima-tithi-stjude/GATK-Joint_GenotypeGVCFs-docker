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
  
  REFERENCE = params.Build_hg38.REFERENCE
  KNOWN_SITE1 = params.Build_hg38.KNOWN_SITE1
  KNOWN_SITE2 = params.Build_hg38.KNOWN_SITE2
  KNOWN_SITE3 = params.Build_hg38.KNOWN_SITE3
  EXOME_BED = params.Build_hg38.EXOME_BED
  INTERVALS = params.Build_hg38.INTERVALS

  hapmap = params.Build_hg38.hapmap 
  omni   = params.Build_hg38.omni
  thousandgenome_snp = params.Build_hg38.thousandgenome_snp
  dbSNP = params.Build_hg38.dbSNP

  mills_indel=params.Build_hg38.mills_indel
  axiom=params.Build_hg38.axiom
  
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
  project          : ${params.project}
  GVCF_filelist    : ${params.GVCF_filelist}
  outdir           : ${params.outdir}
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
  dbSNP            :$dbSNP
  mills_indel      :$mills_indel
  axiom            :$axiom 
  withChr          : $withChr
 """
 // exit 0    
}

def helpMessage() {
  log.info """
        Welocme to run Nextflow Pipeline GenotypeGVCFs.nf 
        Usage:
        A typical command for running the pipeline is as follows:
        nextflow run GenotypeGVCFs.nf -profile cluster_singularity -w ~/Nextflow_work --outdir ./
        A command with more configurable arguments can be           
        nextflow run GenotypeGVCFs.nf -profile cluster_singularity -w ~/Nextflow_work --GVCF_filelist ./nova_GVCF_list.txt --project Ptest --genome_build hg38 --data_type WHOLE_GENOME vqsr_recal Y --outdir ./

        Configurable arguments:
        --GVCF_filelist                  GVCF file list that you want to do the joint cohort calling analysis     
        --project                        Project name/atlas that used to distinguish other analysis           
        --genome_build                   Reference genome version. hg38(default) | hg19 | b37 
        --data_type                      WHOLE_GENOME(default) | EXOME
        --vqsr_recal                     VQSR or not. Y(default) | N
        --ASOption                       Allele Specific Option. AS(Default) | (Empty)
        --popType                        Population datatype, InbreedingCoeff used or not. notPOP(Default) | POP
        --VQSR_indelTranche              VQSR INDEL Tranche level for truth sensitivity filter
        --VQSR_snpTranche                VQSR SNP Tranche level for truth sensitivity filter
        --outdir                         The directory for the pipeline output
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
			   .fromPath(params.GVCF_filelist)               
               .splitText()
               .map{it.replaceFirst(/\n/, "")}
               .map{ file(it) }   //map the file path string into file object, then can extract the file information.  

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
    memory { (32.GB * task.attempt) } // First attempt 32GB, second 64GB, etc
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
    --only-output-calls-starting-in-intervals  \
    --allow-old-rms-mapping-quality-annotation-data  \
    -ploidy 2 
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
    set Chr, file("${params.project}.${Chr}.cohort.vcf.gz"), file("${params.project}.${Chr}.cohort.vcf.gz.tbi") into (cohort_chr_vcf_ch, vqsr_chr_vcf_ch, chr_vcf_recal_ch)

    script:    
    """
    gatk --java-options '-Xmx24g -Xms24g' \
    MergeVcfs \
    ${vcfs.collect{"--INPUT $it " }.join()} \
    -O ${params.project}.${Chr}.cohort.vcf.gz
    """
}

/**************************************************************************************************************************************************
 This process is responsible for preprocess before VQSR recalibration. Each chromosome level joint called vcf file will be preprocessing by getting 
rid of the genotypes, and sorting to generate the index file. The output channel need to fork into three for VQSR recalibration trainning as the type
of snp, and indel respectively, and then the final VQSR applying.  The whole VQSR process can be by passed by params.vqsr_recal      
**************************************************************************************************************************************************/
VQSR_Preprocess_ch =(params.vqsr_recal=="Y")? vqsr_chr_vcf_ch : Channel.empty()
process VQSR_preprocess {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.outdir}/${params.project}/VQSR_preprocess", mode: 'copy', overwrite: true
//  cpus 1
    input:
    set val(Chr), file (chr_vcf), file (chr_vcftbi) from VQSR_Preprocess_ch

    output:
    set Chr, file("${params.project}.${Chr}.vcf.gz"), file("${params.project}.${Chr}.vcf.gz.tbi") into (chr_vcf_indel_ch, chr_vcf_snp_ch)   

    when:
    !(Chr =~ /M/)     //filtered out chrM related
    
    script:
    """
    bcftools view -G -Ou $chr_vcf | bcftools  sort -Oz -o ${params.project}.${Chr}.vcf.gz 
    tabix -p vcf ${params.project}.${Chr}.vcf.gz 
    """
}

/************************************************************************************************************************************************
This process is responsible for recabilicating the short INDEL variants from the chromosomel level vcf files
This process depend on the VQSR_Preprocess and can be by passed.
*************************************************************************************************************************************************/
process VQSR_INDEL_VariantRecalibrator {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.outdir}/${params.project}/VQSR_INDEL_VariantRecalibrator", mode: 'copy', overwrite: true
	
    input:
	set Chrs, file(chrosome_vcfs), file(chrosome_vcftbis) from chr_vcf_indel_ch.groupTuple(by:4)
    set file(genome), file(faidx), file(fadict) from Ref_faidx_vqsr_indel_ch
  
    val dbSNP from dbSNP   
    val DBsnp from KNOWN_SITE2                     
    val mills_indel from mills_indel 
    val axiom from  axiom
    val DPAnnotation //from DPAnnotation
    val ICAnnotation //from ICAnnotation    
    
	output:
    set file("${params.project}.indel.recal"),file("${params.project}.indel.recal.idx"),file("${params.project}.indel.tranches"), file("${params.project}.indel.plots.R.pdf") into indel_recal_ch

    script:
	"""
    gatk --java-options "-Xmx24g -Xms24g" \
    VariantRecalibrator \
    -R ${genome} \
    ${params.ASOption} \
    -resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_indel} \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP} \
    --max-gaussians 4 \
    -an QD \
    ${DPAnnotation} \
    -an FS \
    -an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
    ${ICAnnotation} \
    -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 \
    -tranche 98.8 -tranche 98.6 -tranche 98.4 -tranche 98.2 -tranche 98.0 -tranche 97.8 \
	-tranche 97.6 -tranche 97.4 -tranche 97.2 -tranche 97.0 -tranche 96.8 -tranche 96.6 \
	-tranche 96.4 -tranche 96.2 -tranche 96.0 -tranche 95.8 -tranche 95.6 -tranche 95.4 \
	-tranche 95.2 -tranche 95.0 -tranche 94.8 -tranche 94.6 -tranche 94.4 -tranche 94.2 \
	-tranche 94.0 -tranche 93.8 -tranche 93.6 -tranche 93.4 -tranche 93.2 -tranche 93.0 \
	-tranche 92.8 -tranche 92.6 -tranche 92.4 -tranche 92.2 -tranche 92.0 -tranche 91.8 \
	-tranche 91.6 -tranche 91.4 -tranche 91.2 -tranche 91.0 -tranche 90.8 -tranche 90.6 \
	-tranche 90.4 -tranche 90.2 -tranche 90.0 \
    -mode INDEL \
     ${chrosome_vcfs.collect {"-V $it "}.join()} \
    -O ${params.project}.indel.recal \
    --tranches-file ${params.project}.indel.tranches \
    --rscript-file ${params.project}.indel.plots.R
	"""
}

/************************************************************************************************************************************************
This process is responsible for recabilicating the SNP variants from the output chromosomel level vcf files
This process depend on the VQSR_Preprocess and can be by passed.    
*************************************************************************************************************************************************/
process VQSR_SNP_VariantRecalibrator {
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.outdir}/${params.project}/VQSR_SNP_VariantRecalibrator", mode: 'copy', overwrite: true
 	
    input:
	set Chrs, file(chrosome_vcfs), file(chrosome_vcftbis) from chr_vcf_snp_ch.groupTuple(by:4)
    set file(genome), file(faidx), file(fadict) from Ref_faidx_vqsr_snp_ch
    val hapmap from hapmap
    val omni from omni                  
    val thousandgenome_snp from thousandgenome_snp
    val dbSNP from dbSNP 
    val DPAnnotation // from DPAnnotation
    val MQAnnotation // from MQAnnotation
    val ICAnnotation // from ICAnnotation                   
  
	output:
    set file("${params.project}.snp.recal"),file("${params.project}.snp.recal.idx"),file("${params.project}.snp.tranches"),file("${params.project}.snp.plots.R.pdf"), file("${params.project}.snp.tranches.pdf") into snp_recal_ch

    script:
	"""
    gatk --java-options "-Xmx24g -Xms24g" \
    VariantRecalibrator \
    -R ${genome} \
    ${params.ASOption} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
	-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $thousandgenome_snp \
	-resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $dbSNP \
    --max-gaussians 6 \
    -an QD \
    ${DPAnnotation} \
    -an FS \
    -an SOR \
    ${MQAnnotation} \
    -an MQRankSum \
    -an ReadPosRankSum \
    ${ICAnnotation} \
    -tranche 100.0 -tranche 99.9 -tranche 99.8 -tranche 99.7 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.2 -tranche 99.1 -tranche 99.0 \
    -tranche 98.8 -tranche 98.6 -tranche 98.4 -tranche 98.2 -tranche 98.0 -tranche 97.8 \
	-tranche 97.6 -tranche 97.4 -tranche 97.2 -tranche 97.0 -tranche 96.8 -tranche 96.6 \
	-tranche 96.4 -tranche 96.2 -tranche 96.0 -tranche 95.8 -tranche 95.6 -tranche 95.4 \
	-tranche 95.2 -tranche 95.0 -tranche 94.8 -tranche 94.6 -tranche 94.4 -tranche 94.2 \
	-tranche 94.0 -tranche 93.8 -tranche 93.6 -tranche 93.4 -tranche 93.2 -tranche 93.0 \
	-tranche 92.8 -tranche 92.6 -tranche 92.4 -tranche 92.2 -tranche 92.0 -tranche 91.8 \
	-tranche 91.6 -tranche 91.4 -tranche 91.2 -tranche 91.0 -tranche 90.8 -tranche 90.6 \
	-tranche 90.4 -tranche 90.2 -tranche 90.0 \
    -mode SNP \
    ${chrosome_vcfs.collect {"-V $it "}.join()} \
    -O ${params.project}.snp.recal \
    --tranches-file ${params.project}.snp.tranches \
    --rscript-file ${params.project}.snp.plots.R
	"""
}

/************************************************************************************************************************************************
This process is responsible for applying the recalibration results to the output chromosomel level vcf files
This process depend on the VQSR_INDEL_VariantRecalibrator and VQSR_SNP_VariantRecalibrator, and can be by passed.  
*************************************************************************************************************************************************/
process VQSR_ApplyRecalibration {
   	container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.outdir}/${params.project}/VQSR_ApplyRecalibration", mode: 'copy', overwrite: true

    input:
    set file(genome), file(faidx), file(fadict) from Ref_faidx_vqsr_apply_ch 
    set Chr, file(chr_vcf), file(chr_vcf_tbi) from chr_vcf_recal_ch	   
    set file(indel_recal),file(indel_recal_idx),file(indel_tranches) from indel_recal_ch.first()
    set file(snp_recal),file(snp_recal_idx),file(snp_tranches) from snp_recal_ch.first()

	output:
    set Chr, file("${Chr}.snp.indel.recalibrated.vcf.gz"),file("${Chr}.snp.indel.recalibrated.vcf.gz.tbi") into vqsr_apply_chr_vcf_ch

    when:
    !(Chr =~ /M/)     //filtered out chrM related

    script:
	"""
    gatk --java-options "-Xmx24g -Xms24g" \
    ApplyVQSR \
    -R ${genome} \
    -V ${chr_vcf} \
    -O ${Chr}.indel.recalibrated.vcf.gz \
    ${params.ASOption} \
    --truth-sensitivity-filter-level ${params.VQSR_indelTranche} \
    --recal-file ${indel_recal} \
    --tranches-file ${indel_tranches} \
    -mode INDEL

    gatk --java-options "-Xmx24g -Xms24g" \
    ApplyVQSR \
    -R ${genome} \
    -V ${Chr}.indel.recalibrated.vcf.gz \
    -O ${Chr}.snp.indel.recalibrated.vcf.gz \
    ${params.ASOption} \
    --truth-sensitivity-filter-level ${params.VQSR_snpTranche} \
    --recal-file ${snp_recal} \
    --tranches-file ${snp_tranches} \
    -mode SNP		
	"""
}	

/*********************************************************************************************************************************************************************
  This process is responsible for merging the chrosome level vcf files into a whole vcf files.   
  The input chrosome level vcf files can be from the cohort joint calling result or the VQSR final result.  
*********************************************************************************************************************************************************************/
project_chr_vcf_in =(params.vqsr_recal=="Y")? vqsr_apply_chr_vcf_ch : cohort_chr_vcf_ch
process Merge_Chromosome_To_Project_VCF{
    container 'stithi/gatk_jointgenotype_gvcfs:v1'
    publishDir "${params.outdir}/${params.project}", mode: 'copy', overwrite: true
	  
    input:
    set Chr, file(chrosome_vcfs), file(chrosome_vcftbis) from project_chr_vcf_in.groupTuple(by:4)   // just dump each item into a list just and need to avoid the default key  
    output:
    set file("${params.project}.vcf.gz") , file("${params.project}.vcf.gz.tbi") into project_vcf_ch
    script:
    """
    gatk --java-options "-Xmx24g -Xms24g" \
    MergeVcfs \
    ${chrosome_vcfs.collect{"--INPUT $it " }.join()} \
    -O ${params.project}.vcf.gz 
    """ 
}

