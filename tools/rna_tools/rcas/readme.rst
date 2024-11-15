Important notes
----------------

- the dependencies of the tool are handled by the conda dependency resolver

- to install the dependencies, the galaxy administrator needs to enable the conda channels - bioconda, r, conda-forge, defaults

- by default, the conda installation includes hg19

- to enable analysis for hg38, mm9, mm10, ce10, or dm3, the administrator needs to separately install the relevant packages from the bioconda channel

    - bioconductor-bsgenome.hsapiens.ucsc.hg38
    - bioconductor-bsgenome.mmusculus.ucsc.mm9
    - bioconductor-bsgenome.mmusculus.ucsc.mm10
    - bioconductor-bsgenome.celegans.ucsc.ce10
    - bioconductor-bsgenome.dmelanogaster.ucsc.dm3
    - bioconductor-org.mm.eg.db
    - bioconductor-org.ce.eg.db
    - bioconductor-org.dm.eg.db

