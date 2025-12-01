/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// base
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'
// alignmen
include { STAR_GENOMEGENERATE   } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN            } from '../modules/nf-core/star/align/main'
// quantification
include { SALMON_INDEX          } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT          } from '../modules/nf-core/salmon/quant/main'
// genotyping
include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
// utils
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nf_rna_pipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NF_RNA_PIPELINE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // MODULE: Run FastQC -default is skipped with `skip_tools`
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // PARSE INPUT: Parse samplesheet to separate FASTQ and BAM inputs
    //
    ch_samplesheet
        .branch { meta, fastq1, fastq2, bam ->
            fastq: fastq1 && fastq1 != ""
                return [ meta, fastq2 && fastq2 != "" ? [fastq1, fastq2] : [fastq1] ]
            bam: bam && bam != ""
                return [ meta, bam ]
        }
        .set { ch_input }
    
    //
    // ALIGNMENT (if FQs provided as input)
    //
    // STAR inputs
    ch_star_index = Channel.fromPath(params.star_index).map { [ [:], it ] }
    ch_gtf        = Channel.fromPath(params.gtf).map { [ [:], it ] }
    
    STAR_ALIGN(
        ch_input.fastq,                      // tuple val(meta), path(reads)
        ch_star_index,                       // tuple val(meta2), path(index)
        ch_gtf,                              // tuple val(meta3), path(gtf)
        params.star_ignore_sjdbgtf?: false,  // val star_ignore_sjdbgtf
        params.seq_platform ?: '',           // val seq_platform
        params.seq_center ?: ''              // val seq_center
    )
    
    ch_aligned_bam = STAR_ALIGN.out.bam
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
    
    // Combine BAMs in case fq and bams provided as input
    ch_all_bams = ch_aligned_bam.mix(ch_input.bam)

    //
    // QUANTIFICATION with salmon
    //
    // Create salmon index if needed
    ch_salmon_index = Channel.empty()
    if (params.transcriptome && !params.salmon_index) {
        ch_transcriptome = Channel.fromPath(params.transcriptome)
        SALMON_INDEX(ch_transcriptome)
        ch_salmon_index = SALMON_INDEX.out.index
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())
    } else if (params.salmon_index) {
        ch_salmon_index = Channel.fromPath(params.salmon_index)
    }

    ch_transcript_fasta = params.transcriptome ? 
        Channel.fromPath(params.transcriptome) : Channel.empty()

    if (params.quant_mode.contains('alignment') && ch_all_bams) {
        // Alignment mode with BAM files
        SALMON_QUANT(
            ch_all_bams.map { meta, bam -> [ meta, [bam] ] }, // tuple val(meta), path(reads)
            ch_salmon_index,                                   // path index
            ch_gtf,                                           // path gtf
            ch_transcript_fasta,                              // path transcript_fasta
            true,                                             // val alignment_mode
            'A'                                               // val lib_type
        )
    } else {
        // Mapping mode with FASTQ files
        SALMON_QUANT(
            ch_input.fastq,              // tuple val(meta), path(reads)
            ch_salmon_index,             // path index
            ch_gtf,                      // path gtf
            ch_transcript_fasta,         // path transcript_fasta
            false,                       // val alignment_mode
            'A'                          // val lib_type
        )
    }
    ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())


    //
    // GENOTYPING with GATK's HaplotypeCaller
    //
    // Prepare reference channels
    ch_fasta = Channel.fromPath(params.fasta).map { [ [:], it ] }
    ch_fai   = Channel.fromPath(params.fasta_fai).map { [ [:], it ] }
    ch_dict  = Channel.fromPath(params.dict).map { [ [:], it ] }

    // Prepare input for GATK4_HAPLOTYPECALLER
    ch_bam_for_hc = ch_all_bams.map { meta, bam ->
        // nf-core module expects: tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
        [ meta, bam, [], [], [] ]
    }

    GATK4_HAPLOTYPECALLER(
        ch_bam_for_hc,                           // tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
        ch_fasta,                                // tuple val(meta2), path(fasta)
        ch_fai,                                  // tuple val(meta3), path(fai)
        ch_dict,                                 // tuple val(meta4), path(dict)
        Channel.empty().map { [ [:], [] ] },     // tuple val(meta5), path(dbsnp)
        Channel.empty().map { [ [:], [] ] }      // tuple val(meta6), path(dbsnp_tbi)
    )
    
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())
    
    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf_rna_pipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
