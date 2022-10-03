#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/genomeannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/genomeannotator
    Website: https://nf-co.re/genomeannotator
    Slack  : https://nfcore.slack.com/channels/genomeannotator
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GeneidX PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// nextflow run main.nf -profile docker --assembly data/SampleGenome.fa --outdir output/ --rm_lib data/GCA_000180735.1.repeatmodeler.fa --prot_file output/proteins/UniRef90.6658.3+.fa.gz

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

OutputFolderProteinDBs = "${OutputFolder}/proteins"
OutputFolderSpecies = "${OutputFolder}/species"
OutputFolderSpeciesTaxid = "${OutputFolder}/species/${params.taxid}"
// fasta.gz
// fasta.gz.gzi
// fasta.gz.fai
// gff3.gz
// gff3.gz.tbi
// hsp.gff
// param

// paramOutputFolder = "${params.output}/params"

// genoom = file(params.genome)


/*
 * Defining the module / subworkflow path, and include the elements
 */
subwork_folder = "${projectDir}/subworkflows"

include { UncompressFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolder)

include { prot_down_workflow } from "${subwork_folder}/getProteins" addParams(OUTPUT: OutputFolderProteinDBs,
 LABEL:'singlecpu')

include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolderProteinDBs,
 LABEL:'fourcpus')

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'fourcpus')

include { matchAssessment } from "${subwork_folder}/getTrainingSeq" addParams(OUTPUT: OutputFolder,
 LABEL:'singlecpu')

include { param_selection_workflow } from "${subwork_folder}/getParams" addParams(OUTPUT: OutputFolder,
 LABEL:'singlecpu')

include { param_value_selection_workflow } from "${subwork_folder}/getParamsValues" addParams(OUTPUT: OutputFolder,
 LABEL:'singlecpu')

include { creatingParamFile } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

include { creatingParamFile_frommap } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams( LABEL:'singlecpu' )

include { prep_concat } from "${subwork_folder}/prepare_concatenation" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

include { concatenate_Outputs_once } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

include { gff3addInfo } from "${subwork_folder}/addMatchInfo" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

// compress and index fastas to be stored and published to the cluster
include { compress_n_indexFASTA } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid,
 LABEL:'singlecpu')

// compress and index gff3s to be stored and published to the cluster
include { gff34portal } from "${subwork_folder}/tools" addParams(OUTPUT: OutputFolderSpeciesTaxid)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEANNOTATOR } from './workflows/genomeannotator'

//
// WORKFLOW: Run main nf-core/genomeannotator analysis pipeline
//
// workflow NFCORE_GENOMEANNOTATOR {
//     GENOMEANNOTATOR ()
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {

    parameter_location = file(params.parameter_path)

    uncompressed_genome = GENOMEANNOTATOR()


    // none of the returned objects is used by downsteam processes
    compress_n_indexFASTA(uncompressed_genome)


    // if proteins_file provided use proteins file
    // else, use taxon to download the proteins_file
    // both conditions are evaluated inside the execution of this workflow
    if (params.prot_file) {
      proteins_file = file(params.prot_file)
    } else {
      proteins_file = prot_down_workflow(params.taxid,
                                         params.proteins_lower_lim,
                                         params.proteins_upper_lim)
    }

    // Build protein database for DIAMOND
    protDB = build_protein_DB(proteins_file)



    // Run DIAMOND to find matches between genome and proteins
    hsp_found = alignGenome_Proteins(protDB, uncompressed_genome)


    // Automatic computation of the parameter file
    new_mats = matchAssessment(uncompressed_genome, hsp_found,
                                    params.match_score_min,
                                    params.match_ORF_min,
                                    params.intron_margin,
                                    params.min_intron_size,
                                    params.max_intron_size)


    // if sites matrices provided, use them
    // else, use taxon to get the closest geneid param file
    if (params.acceptor_pwm) {
      acc_pwm = params.acceptor_pwm
      don_pwm = params.donor_pwm
      sta_pwm = params.start_pwm
      sto_pwm = params.stop_pwm
    } else {
      param_file_sel = param_selection_workflow(params.taxid, 0, parameter_location)
      acc_pwm = param_file_sel.acceptor_pwm
      don_pwm = param_file_sel.donor_pwm
      sta_pwm = param_file_sel.start_pwm
      sto_pwm = param_file_sel.stop_pwm
    }

    para_vals = param_value_selection_workflow(params.taxid, 0,
                                              parameter_location,
                                              params.maps_param_values)

    new_param = creatingParamFile_frommap(
                                          params.taxid,

                                          para_vals.params_map,

                                          sta_pwm,
                                          acc_pwm,
                                          don_pwm,
                                          sto_pwm,

                                          new_mats.ini_comb,
                                          new_mats.trans_comb,

                                          params.general_gene_params,

                                          hsp_found

                                          )


    // Run Geneid
    predictions = geneid_WORKFLOW(uncompressed_genome,
                                  new_param,
                                  hsp_found)

    // Prepare concatenation
    // This will initialize the final GFF3 file
    output_file = prep_concat(proteins_file, uncompressed_genome)



    // Run concatenation of individual GFF3 files
    final_output = concatenate_Outputs_once(predictions.collect(), output_file)

    // Add information about the proteins to the final GFF3
    labelled_output = gff3addInfo(final_output, hsp_found)


    // If proteins from UniRef, add GO terms to the GFF3
    if (params.source_uniprot){
      // func_labelled_output = addGOterms(labelled_output)
      // gff34portal(func_labelled_output)
      gff34portal(labelled_output)
    } else {
      // fix gff3 file and compress it for the portal
      gff34portal(labelled_output)
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
