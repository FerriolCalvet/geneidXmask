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

/*
 * Defining the output folders.
 */
OutputFolder = "${params.output}"

OutputFolderInternal = "${OutputFolder}/internal"
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

include { prot_down_workflow } from "${subwork_folder}/getProteins" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { filter_Fasta_by_length } from "${subwork_folder}/filter_fasta" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { build_protein_DB } from "${subwork_folder}/build_dmnd_db" addParams(OUTPUT: OutputFolderInternal,
  LABEL:'fourcpus')

include { alignGenome_Proteins } from "${subwork_folder}/runDMND_BLASTx" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'fourcpus')

include { geneid_WORKFLOW } from "${subwork_folder}/geneid" addParams( LABEL:'singlecpu' )

include { prep_concat } from "${subwork_folder}/prepare_concatenation" addParams(OUTPUT: OutputFolderSpeciesTaxid,
  LABEL:'singlecpu')

include { concatenate_Outputs } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { concatenate_Outputs_once } from "${subwork_folder}/geneid_concatenate" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { matchAssessment } from "${subwork_folder}/getTrainingSeq" addParams(OUTPUT: OutputFolder,
  LABEL:'singlecpu')

include { creatingParamFile } from "${subwork_folder}/modifyParamFile" addParams(OUTPUT: OutputFolderSpeciesTaxid,
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


    // if proteins_file provided use proteins file
    // else, use taxon to download the proteins_file
    // both conditions are evaluated inside the execution of this workflow
    if (params.prot_file) {
      proteins_file = file(params.prot_file)
    } else {
      proteins_file = prot_down_workflow(params.taxid)
    }

    // Build protein database for DIAMOND
    protDB = build_protein_DB(proteins_file)

    // // Run DIAMOND to find matches between genome and proteins
    // rep_matched_found = alignGenome_Proteins(protDB, repeats_data)


    uncompressed_genome = GENOMEANNOTATOR()


    // none of the returned objects is used by downsteam processes
    compress_n_indexFASTA(uncompressed_genome)

    // Remove contigs that are too small (user chooses threshold)
    // filtered_genome = filter_Fasta_by_length(uncompressed_genome, params.min_seq_length)





    // // Build protein database for DIAMOND
    // protDB = build_protein_DB(proteins_file)


    // Run DIAMOND to find matches between genome and proteins
    hsp_found = alignGenome_Proteins(protDB, uncompressed_genome)


    // Automatic computation of the parameter file
    new_mats = matchAssessment(uncompressed_genome, hsp_found,
                                    params.match_score_min,
                                    params.match_ORF_min,
                                    params.intron_margin,
                                    params.min_intron_size,
                                    params.max_intron_size)


    new_param = creatingParamFile(uncompressed_genome,
                                  params.no_score,
                                  params.site_factor,
                                  params.exon_factor,
                                  params.hsp_factor,
                                  params.exon_weight,

                                  params.min_intron_size_geneid,
                                  params.max_intron_size_geneid,

                                  params.start_pwm,
                                  params.acceptor_pwm,
                                  params.donor_pwm,
                                  params.stop_pwm,

                                  new_mats.ini_comb,
                                  new_mats.trans_comb
                                  )



    // Run Geneid
    predictions = geneid_WORKFLOW(uncompressed_genome, new_param, hsp_found)

    // Prepare concatenation
    // This will initialize the final GFF3 file
    output_file = prep_concat(proteins_file, uncompressed_genome)



    // Run concatenation of individual GFF3 files
    // final_output = concatenate_Outputs(predictions, output_file)
    final_output = concatenate_Outputs_once(predictions.collect(), output_file)


    // fix gff3 file and compress it for the portal
    gff34portal(final_output.last())

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/