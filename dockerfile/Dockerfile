# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Install required Bioconductor package
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'install.packages("argparse")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("ggrepel")'
