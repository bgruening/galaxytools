Important notes
----------------

- the dependencies of the tool are handled by the conda dependency resolver

- to install the dependencies, the galaxy administrator needs to configure to ensure the conda channels - bioconda, r, conda-forge, defaults

- by default, the conda installation includes hg19

- to enable analysis for mm9, ce10, dm3, the administrator needs to separately install the relevant packages from the channel bioconda
    - bioconductor-bsgenome.dmelanogaster.ucsc.dm3
    - bioconductor-bsgenome.celegans.ucsc.ce10
    - bioconductor-bsgenome.mmusculus.ucsc.mm9
    - bioconductor-org.mm.eg.db
    - bioconductor-org.ce.eg.db
    - bioconductor-org.dm.eg.db
