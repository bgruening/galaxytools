---
title: Dependency Handling
layout: default
---

|

Many tools using common libraries to base their work on.
For example bowtie2, numpy, matplotlib or GraphicsMagic are widly used from several projects. 
The Galaxy Toolshed supports handling of dependencies and if your programm uses some external library that meight be of interest by other developers, you should create an `orphan tool dependency`_.
Orphan tool dependencies are packages that do not contain any tools, but definitions for libraries or binaries that can be easily used by tool-developers. 
You can find several examples in the ChemicalToolBoX_. Following a list of `orphan tool dependency`_ we have developed during the time.


- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_atlas_3_10 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_boost_1_53
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_bzlib_1_0
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_chemfp_1_1
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_eigen_2_0
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_eigen_3_1
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_freetype_2_4
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_graphicsmagick_1_3 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_lapack_3_4 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_matplotlib_1_2
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_numpy_1_7 
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_openbabel_2_3
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_rdkit_2012_12
- http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_rdkit_2013_03


.. _`orphan tool dependency`: http://wiki.galaxyproject.org/DefiningRepositoryDependencies
.. _ChemicalToolBoX: /galaxytools/projects/chemicaltoolbox/
