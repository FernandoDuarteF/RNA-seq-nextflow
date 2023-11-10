## RNA-seq

### Processing

This pipeline is for paired end reads RNA-seq sequences. The final output is matrix called ``raw_counts_matrix.txt`` with counts for each sample, located in ``$projectDir/results/folder``.

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

In ``nextflow.config`` is not necessary to indicate a default container if it's not necessary. If a linux environment is needed for a specific process without a container, use ``biocorecrg/debian-perlbrew:latest`` from quay.io.

For trimmomatic process in ``main.nf``, choose extraparameters for adapter and quality trimming. ``ILLUMINACLIP:${adaptors}:2:30:10`` is placed as an example.

Final output is ``*ReadsPerGene.out.tab`` count files.

If strandness is known, the MATRIX process can be used for building the count matrix from the STAR outpt. Use the correct column number (2, 3, or 4, according to strandness) from the ``*ReadsPerGene.out.tab`` files (refer to STAR output explanation) as an input parameter. Refer to STAR manual if unclear.

### Downstream analysis

You can use the dockerfile to build the docker image. After bulding it, you can run the DESeq2 script in a docker cointainer:

```
# Build docker image
docker build -t docker_image:tag
# Mount folder with script and input files into a new container
docker run --detach --volume $(pwd):/scratch --name deseq2_container docker_image:tag tail -f /dev/null
# Run script inside container
docker run -it deseq2_container /bin/bash
Rscript deseqV2.R -t STAR -g gene
```

Do ``Rscript deseqV2.R --help`` for argument options. The script above uses default arguments.

Bear in mind that this script was made with single factor in mind (Genotype), and two levels (wt and ko, view example metadata example). Make any necessary changes to adjust the script to your experimental design.

**IN PROGRESS**

