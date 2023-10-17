###RNA-seq

##Processing

This pipeline is for paired end reads RNA-seq sequences. The final output is matrix called ``raw_counts_matrix.txt`` with counts for each sample, located in ``$projectDir/results/folder``.

Run pipeline using:

```
nextflow run main.nf -profile crg -resume -bg > log
```

Indicate paths to fastq and index files using ``params.config``, as well as results folder. Modify cluster options in ``nextflow.config`` according to the cluster executor. You can change ``params.config`` container to whatever fits best. This container is only used to run a bash script.

In ``nextflow.config`` is not necessary to indicate a default container if it's not necessary. If a linux environment is needed for a specific process without a container, use ``biocorecrg/debian-perlbrew:latest`` from quay.io. 
