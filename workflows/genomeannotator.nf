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


OutputFolder = "${params.output}"

/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows"

include { repeat_down_workflow } from "${subwork_folder}/getRepeats" addParams(OUTPUT: OutputFolder,
 LABEL:'singlecpu')



// Set relevant input channels
if (params.rm_lib) {
  ch_repeats_status = 1
  ch_repeats = Channel.fromPath(file(params.rm_lib, checkIfExists: true))
} else {
  ch_repeats_status = 0
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
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

  if (!ch_repeats_status) {
    data_location = file(params.repeats_data_path)
    ch_repeats = repeat_down_workflow(params.taxid,
                                        data_location)
  }


  //
  // SUBWORKFLOW: Validate and pre-process the assembly
  //
  ASSEMBLY_PREPROCESS(
      ch_genome
  )
  ch_versions = ch_versions.mix(ASSEMBLY_PREPROCESS.out.versions)


  //
  // SUBWORKFLOW: Repeat modelling if this option is requested
  //
  if (params.repeat_modeler) {
     REPEATMODELER(
        ASSEMBLY_PREPROCESS.out.fasta
     )
     ch_repeats = REPEATMODELER.out.fasta.map {m,fasta -> fasta}

     REPEATMASKER(
        ASSEMBLY_PREPROCESS.out.fasta,
        ch_repeats,
        false,
        ch_rm_db
     )
     ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
     ch_genome_rm = REPEATMASKER.out.fasta

  } else {
    // otherwise directly mask using the selected repeats

    REPEATMASKER(
       ASSEMBLY_PREPROCESS.out.fasta,
       ch_repeats,
       false,
       ch_rm_db
    )
    ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
    ch_genome_rm = REPEATMASKER.out.fasta

  }

  // //
  // // MODULE: Repeatmask the genome; if a repeat species is provided, use it
  // //         - else the repeats in FASTA format
  // // NOT REALLY SURE HOW THE REPEAT SPECIES THING WORKS...
  // if (params.rm_species) {
  //    REPEATMASKER(
  //       ASSEMBLY_PREPROCESS.out.fasta,
  //       ch_repeats,
  //       params.rm_species,
  //       ch_rm_db
  //    )
  //    ch_versions = ch_versions.mix(REPEATMASKER.out.versions)
  //    ch_genome_rm = REPEATMASKER.out.fasta
  // }


  emit:
  ch_genome_rm
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
