/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGenomeannotator.initialise(params, log)

def checkPathParamList = [ params.assembly ]
for (param in checkPathParamList) {
  if (param) {
    file(param, checkIfExists: true)
  }
}

// Check mandatory parameters
if (params.assembly) {
  ch_genome = file(params.assembly, checkIfExists: true)
} else {
  exit 1, 'No assembly specified!'
}

// Set relevant input channels
if (params.rm_lib) {
  ch_repeats = Channel.fromPath(file(params.rm_lib, checkIfExists: true))
} else {
  ch_repeats = Channel.fromPath("${workflow.projectDir}/assets/repeatmasker/repeats.fa")
}

if (params.rm_db) {
  ch_rm_db = file(params.rm_db)
} else {
  ch_rm_db = Channel.empty()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// ch_rfam_cm = file("${workflow.projectDir}/assets/rfam/14.2/Rfam.cm.gz", checkIfExists: true)
// ch_rfam_family = file("${workflow.projectDir}/assets/rfam/14.2/family.txt.gz", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { ASSEMBLY_PREPROCESS } from '../subworkflows/local/assembly_preprocess'
include { REPEATMASKER } from '../subworkflows/local/repeatmasker'
// include { BUSCO_QC } from '../subworkflows/local/busco_qc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
// include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
// include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { REPEATMODELER } from '../modules/local/repeatmodeler'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GENOMEANNOTATOR {

    ch_versions = Channel.empty()
    ch_repeats_lib = Channel.empty()
    ch_genome_rm = Channel.empty()


    //
    // SUBWORKFLOW: Validate and pre-process the assembly
    //
    ASSEMBLY_PREPROCESS(
        ch_genome
    )
    ch_versions = ch_versions.mix(ASSEMBLY_PREPROCESS.out.versions)


    //
    // SUBWORKFLOW: Repeat modelling if no repeats are provided
    //
    if (!params.rm_lib && !params.rm_species) {
       REPEATMODELER(
          ASSEMBLY_PREPROCESS.out.fasta
       )
       ch_repeats = REPEATMODELER.out.fasta.map {m,fasta -> fasta}
    }

    //
    // MODULE: Repeatmask the genome; if a repeat species is provided, use that - else the repeats in FASTA format
    if (params.rm_species) {
       REPEATMASKER(
          ASSEMBLY_PREPROCESS.out.fasta,
          ch_repeats,
          params.rm_species,
          ch_rm_db
       )
       ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
       ch_genome_rm = REPEATMASKER.out.fasta
    } else {
       REPEATMASKER(
          ASSEMBLY_PREPROCESS.out.fasta,
          ch_repeats,
          false,
          ch_rm_db
       )
       ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
       ch_genome_rm = REPEATMASKER.out.fasta
    }
    emit:
    ch_genome_rm
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    // if (params.email || params.email_on_fail) {
    //     NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    // }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
