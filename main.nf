include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run ebi-metagenomics/miniqc-pipeline --samplesheet input_file.csv --ref_genome bwa-mem2-indexed.fasta")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Create a new channel of metadata from a sample sheet passed to the pipeline through the --input parameter
ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_input.json"))

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.fromPath("$projectDir/assets/mgnify_logo.png")

/*
~~~~~~~~~~~~~~~~~~
    DBs
~~~~~~~~~~~~~~~~~~
*/
ref_genome         = file(params.ref_genome)
ref_genome_index   = file("${ref_genome.parent}/*.fa.*")

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_BEFORE     } from './modules/fastqc'
include { FASTQC as FASTQC_AFTER      } from './modules/fastqc'
include { MULTIQC                     } from './modules/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/dumpsoftwareversions'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     Subworkflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { QC_AND_MERGE_READS   } from './subworkflows/qc_and_merge_reads'
include { DECONTAMINATION      } from './subworkflows/decontamination'

workflow {

    ch_versions    = Channel.empty()

    ref_genome         = file(params.ref_genome)
    ref_genome_index   = file("${ref_genome.parent}/*.${ref_genome.extension}.*")

    /* Adjust the input structure */
    groupReads = { meta, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta + [single_end: true], [fq1])
        }
        else {
            return tuple(meta + [single_end: false], [fq1, fq2])
        }
    }
    ch_input_runs = ch_input.map(groupReads) // [ meta, [raw_reads] ]

    // --- check input reads quality --- //
    FASTQC_BEFORE( ch_input_runs )

    ch_versions = ch_versions.mix(FASTQC_BEFORE.out.versions)

    // --- quality control ---- //
    QC_AND_MERGE_READS( ch_input_runs )

    ch_versions = ch_versions.mix( QC_AND_MERGE_READS.out.versions )

    // --- decontamination ---- //
    // We need a tuple as the alignment and decontamination module needs the input like that
    DECONTAMINATION( QC_AND_MERGE_READS.out.reads, ref_genome, ref_genome_index )

    ch_versions = ch_versions.mix( DECONTAMINATION.out.versions )

    // --- check filtered reads quality --- //
    FASTQC_AFTER ( DECONTAMINATION.out.decontaminated_reads )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_BEFORE.out.zip.collect{it[1]}.ifEmpty([]) )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_AFTER.out.zip.collect{it[1]}.ifEmpty([]) )
    ch_multiqc_files = ch_multiqc_files.mix( QC_AND_MERGE_READS.out.mqc.map { map, json -> json }.collect().ifEmpty([]) )
    ch_multiqc_files = ch_multiqc_files.mix( CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect() )

    MULTIQC(
         ch_multiqc_files.collect(),
         ch_multiqc_config.toList(),
         ch_multiqc_custom_config.toList(),
         ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

}