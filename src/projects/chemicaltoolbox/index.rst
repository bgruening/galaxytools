---
title: ChemicalToolBoX
layout: default
---


- Installation_
- Tools_

.. _Installation: /galaxytools/projects/chemicaltoolbox/installation
.. _Tools: /galaxytools/projects/chemicaltoolbox/tools


===========================
What's the ChemicalToolBoX?
===========================

The ChemicalToolBoX is a set of tools integrated into the Galaxy-workflow-management system to enable researchers easy-to-use, reproducible, and transparent access to 
cheminformatics libraries and drug discovery tools. It includes standard applications for similarity and 
substructure searches, clustering of compounds, prediction of properties and descriptors, filtering, and many 
other tools that range from drug-likeness classification to fragmentation and fragment-merging.
By combinating the various tools many more powerful applications can be designed.

ChemicalToolBox is based on open-source software, web-accessible, freely available, and easily expandable. 
It can be downloaded and easily deployed locally or on a large scale cluster.

Galaxy
======

`Galaxy <http://galaxyproject.org/>`_ is an open, web-based platform for data intensive research.
All tools can be combined in workflows without any need of programming skills. 
Furthermore the platform can be extended with more tools at any time.
Each tool has its own information about what it does and how the input is supposed to look like.
You can make data available for Galaxy by uploading local files or downloading online content.
Inputfiles, workflowsteps and results are stored in a history where you can view them or reaccess them later.
It is possible to share workflows and histories with other users or make the public available.
Saved workflows can be used with new input files or just to rerun an analyses which ensures repeatability.

Parallelisation
===============

ChemicalToolBoX is capable of accelerating the computation time of resource-intensive processes. 
Large molecule files are automatically splitted in smaller chunks. Each chunk then is processed on a separate core and merged afterwards. 
Everything stays the same for the user but computation time decreases.

As the parallelisation is scalable the job will run on a predefined number of cores. The more cores the faster the processing.

===================
Supported Filetypes
===================

- InChI_
	International Chemical Identifier - developed by the IUPAC_. Representation of a chemical molecule as a string which can include information about the bond, tautomerism, isotope, charge and stereochemistry. Strings are generated following the InChI-algorithm.
- MOL_ & MOL2_
	Single chemical molecule. See SDF.
- SDF_
	Structure-data-file consisting of many MOL-files. Molecules are separated by four Dollar signs ($$$$). Allows the storing metainformation like molecular mass or uniqueID. Developed by MDL Information System (Accelrys_).
- SMILES_
	A line notation using ASCII strings to represent chemical molecules. Information about the charge, isotope or radical can be included besides the stereo information (CIP convention) and the normal bonds. The Simplified Molecular Input Line Entry Specification was developed by Daylight_ Chemical Information System Incorporation.

.. _InChI: http://www.iupac.org/home/publications/e-resources/inchi.html
.. _IUPAC: http://www.iupac.org

.. _MOL: http://en.wikipedia.org/wiki/Chemical_table_file
.. _MOL2: http://openbabel.org/wiki/Sybyl_mol2
.. _SDF: http://accelrys.com/products/informatics/cheminformatics/ctfile-formats/no-fee.php
.. _Accelrys: http://accelrys.com

.. _SMILES: http://daylight.com/smiles/index.html
.. _Daylight: http://daylight.com

All filetypes are interchangable due to three easy converting options:


- the built-in conversion via the pencil icon |pencilicon|

	.. |pencilicon| image:: http://github.com/sbleher/galaxytools/raw/master/chemicaltoolbox/readme/pencilicon.png
	.. image:: http://github.com/sbleher/galaxytools/raw/master/chemicaltoolbox/readme/convert_pencil.jpg

- the Compound Converter tool described in the Tools section

- the import conversion each tool automatically offers

	.. image:: http://github.com/sbleher/galaxytools/raw/master/chemicaltoolbox/readme/internal_conversion.png




Bug Tracker
===========
Have a bug or a feature request? `Please write a new card`_. Before writing a new card, please search for existing issues.

.. _Please write a new card: https://trello.com/b/t9Wr8lSY


Contributing
============
We encourage you to contribute to ChemicalToolBoX! Check out our `Trello board`_ or contact us via e-mail_.

.. _Trello board: https://trello.com/b/t9Wr8lSY
.. _e-mail: (bjoern.gruening@gmail.com)|: 
