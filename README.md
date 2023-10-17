### RNA-seq

## Processing

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

In ``nextflow.config`` is not necessary to indicate a default container if it's not necessary. If a linux environment is needed for a specific process without a container, use ``biocorecrg/debian-perlbrew:latest`` from quay.io.

For trimmomatic process in ``main.nf``, choose extraparameters for adapter and quality trimming. ``ILLUMINACLIP:${adaptors}:2:30:10`` is placed as an example.

Final output is ``*ReadsPerGene.out.tab`` count files.
