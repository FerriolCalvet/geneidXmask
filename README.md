# GeneidX with masking step

In this repository we provide a pipeline for the annotation of any genome starting from any raw genome assembly.

The two main steps are an initial genome masking step that was adapted from the nf-core/genomeannotator pipeline, and then a gene prediction step using GeneidX.

Check our repository of [GeneidX](https://github.com/guigolab/geneidx) for more information on the method and a benchmarking with several species.

Stay tuned for an article with detailed descriptions and feel free to [contact us](mailto:ferriol.calvet@crg.eu) if you are trying it and find any problem.


## Running GeneidX with masking
Having defined the parameters and ensuring that Nextflow and a container technology.

This pipeline requires the uncompressed genome assembly and the taxid of the species to annotate.

`nextflow run ferriolcalvet/geneidXmask -profile <docker/singularity>
                                        --assembly <GENOME>.fa
                                        --taxid <TAXID>
                                        --outdir <OUTPUT_directory>`

or alternatively, clone the repository and then run it (highly recommended)
`git clone https://github.com/FerriolCalvet/geneidXmask.git`
`cd geneidXmask`
`nextflow run main.nf -profile <docker/singularity>
                      --assembly <GENOME>.fa
                      --taxid <TAXID>
                      --outdir <OUTPUT_directory>`


## Before running GeneidXmask
1. **Make sure ( Docker or Singularity ) and Nextflow are installed in your computer.**
  - [Docker](https://docs.docker.com/engine/install/)
  - [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html#)
  - [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Define the input parameters:
  - Uncompressed FASTA file with `.fa` termination.
  - Taxid of the species to annotate or from the closest relative with a described taxid.
  - A repeats file can be optionally provided to the program, otherwise these repeats will be automatically chosen.
  - Also a set of proteins can be provided or otherwise they are automatically chosen.
  - See [GeneidX](https://github.com/guigolab/geneidx) repository for more parameters that can be tuned.


## Work in progress
  - It is recommended to clone the repository and then run the pipeline from there.
  - So far the output of the predictions is not stored in the path indicated when running the pipeline, but in output/species/{taxid of the species}.
  - If you are running the pipeline multiple times, it is recommended that you define a directory for downloading the docker/singularity images to avoid having to download them multiple times. See `singularity.cacheDir variable` in `nextflow.config`.



Follow us on Twitter ([@GuigoLab](https://twitter.com/GuigoLab)) for updates in the article describing our method.
