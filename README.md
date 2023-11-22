## RNA-seq

### Processing

This pipeline is for paired end reads RNA-seq sequences. The final output is matrix called ``raw_counts_matrix.txt`` with counts for each sample, located in ``$projectDir/results/folder``. The pipeline is made to run on a specific cluster using a SGE based executor. Modify profile options in ``nextflow.config`` accordingly.

Run pipeline using:

```
nextflow run main.nf -profile crg -resume -bg > log
```

Indicate paths to fastq and index files using ``params.config``, as well as results folder. Modify cluster options in ``nextflow.config`` according to the cluster executor. You can change ``params.config`` container to whatever fits best. This container is only used to run a bash script or perform other Linux tasks.

Parameters needed:

- fastq.gz files (forward and reverse)
- reference.fasta file for building the index
- annotation.gtf file for building the index
- overhang value. View STAR docummnetation for correct value. FastQC on reads might be need to be ran first.
- adaptors.fasta file (leave empty ("") if adaptor trimming is not necessary)

Note that reference files (fasta and gtf) need to be uncompressed for building the index (for now).

Is not necessary to indicate a default container in ``nextflow.config``, but you can do so in order to run a process that has no containers asigned (not the case here). If a linux environment is needed for a specific process with no containers asigned, use ``biocorecrg/debian-perlbrew:latest`` from quay.io.

For trimmomatic process in ``main.nf``, choose extraparameters for adapter and quality trimming. ``ILLUMINACLIP:${adaptors}:2:30:10`` is placed as an example.

Final output is ``*ReadsPerGene.out.tab`` count files.

If strandness is known, the MATRIX process can be used for building the count matrix from the STAR outpt. Use the correct column number (2, 3, or 4, according to strandness) from the ``*ReadsPerGene.out.tab`` files as an input parameter. Refer to STAR manual if unclear.

**PIPELINE IS STILL IN PROGRESS** Salmon module is to be included. Also, automatically detect read length based on mode using fastq output.

### Downstream analysis

You can use the dockerfile to build the docker image. After bulding it, you can run the DESeq2 script in a docker cointainer:

```
# Build docker image
docker build -t docker_image:tag
# Mount folder with script and input files into a new container
docker run --detach --volume $(pwd):/scratch --name deseq2_container docker_image:tag tail -f /dev/null
# Run script inside container
docker run -it deseq2_container /bin/bash
Rscript deseqV2.R --help
```

Alternatively, you can use singularity to run the script. This might be necessary if you are running the script from the cluster:

```
singularity pull deseq2.sif docker://dernandod95/deseq2:deseqV2
singulairty exec -e deseq2.sif Rscript deseqV2.R --help
```

Run ``Rscript deseqV2.R --help`` for argument options:

```
usage: dockerV2.R [-h] [-st STAR_MATRIX] [-sl SALMON_FOLDER] [-t2g TX2GENE]
                  [-m METADATA] [-t DATA_TYPE] [-g GENE_COUNTS]

options:
  -h, --help            show this help message and exit
  -st STAR_MATRIX, --star_matrix STAR_MATRIX
                        Path to STAR raw counts table
  -sl SALMON_FOLDER, --salmon_folder SALMON_FOLDER
                        Path to salmon count folders
  -t2g TX2GENE, --tx2gene TX2GENE
                        Path to tx2gene table
  -m METADATA, --metadata METADATA
                        Path to metadata table
  -t DATA_TYPE, --data_type DATA_TYPE
                        Run the script on STAR or salmon output [STAR,salmon]
  -g GENE_COUNTS, --gene_counts GENE_COUNTS
                        Output gene counts for specific gene based on tx2gene
                        table gene names (optional)
```

Usage example:

```
# For STAR
singulairty exec -e deseq2.sif Rscript deseqV2.R -st ./raw_counts_matrix.txt -t2g ./tx2gene.gencode.csv -m ./meta_sorted.tsv -t STAR -g Arntl
# For salmon
singulairty exec -e deseq2.sif Rscript deseqV2.R -st ./raw_counts_matrix.txt -t2g ./tx2gene.gencode.csv -m ./meta_sorted.tsv -t STAR -g Arntl
```

Take into account that this script was made with single factor in mind (Genotype), and two levels (wt and ko). View ``DESeq2`` folder for **tx2gene** and **metadata** for examples and format. Make any necessary changes to adjust the script to your experimental design.

**IN PROGRESS. SALMON NEEDS TO BE INCLUDED. ALSO ADD OVEHANG CALCULATION**

* Add RSeQC for gene body coverage for STAR bam output
* Add line to join counts from star output
