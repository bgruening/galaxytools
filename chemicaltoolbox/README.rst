***************
ChemicalToolBoX
***************

Contents
========
- `What's the ChemicalToolBoX?`_
	- `About Galaxy`_
	- Parallelisation_
	- `Supported Filetypes`_

- `Getting Started`_
	- `Learn Galaxy`_
	- Installation_
	- `Admin Account`_
- Toolshed_
- Tools_
- `Bug Tracker`_
- Contributing_

____________________________


.. _Learn Galaxy: http://wiki.galaxyproject.org/Learn
.. _What's the ChemicalToolBoX?

What's the ChemicalToolBoX?
===========================

The ChemicalToolBoX is a set of tools integrated into the Galaxy-workflow-management system to enable researchers easy-to-use, reproducible, and transparent access to 
cheminformatics libraries and drug discovery tools. It includes standard applications for similarity and 
substructure searches, clustering of compounds, prediction of properties and descriptors, filtering, and many 
other tools that range from drug-likeness classification to fragmentation and fragment-merging.
By combinating the various tools many more powerful applications can be designed.

ChemicalToolBox is based on open-source software, web-accessible, freely available, and easily expandable. 
It can be downloaded and easily deployed locally or on a large scale cluster.

.. _About Galaxy:
======
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

.. _Parallelisation:
===============
Parallelisation
===============

ChemicalToolBoX is capable of accelerating the computation time of resource-intensive processes.
Large molecule files are automatically splitted in smaller chunks.
Each chunk then is processed on a separate core and merged afterwards.
Everything stays the same for the user but computation time decreases.

As the parallelisation is scalable the job will run on a predefined number of cores.
The more cores the faster the processing.

.. _Supported Filetypes:
===================
Supported Filetypes
===================

- InChI_
	International Chemical Identifier - developed by the IUPAC_. Representation of a chemical molecule as a string which can include information about the bond, tautomerism, isotope, charge and stereochemistry. Strings are generated following the InChI-algorithm.
- MOL2_
	 A Tripos_ Mol2 file can store a complete representation of a SYBYL molecule.
- MOL_ & SDF_
	Structure-data-file consisting of many MOL-files. Molecules are separated by four Dollar signs ($$$$). Allows the storing metainformation like molecular mass or uniqueID. Developed by MDL Information System (Accelrys_).
- SMILES_
	A line notation using ASCII strings to represent chemical molecules. Information about the charge, isotope or radical can be included besides the stereo information (CIP convention) and the normal bonds. The Simplified Molecular Input Line Entry Specification was developed by Daylight_ Chemical Information System Incorporation.
- and others:
	Special filetypes like the `Open Babel`_ Fastsearch_ index or the Pharmacophore type from `silicos-it` are also supported.

.. _InChI: http://www.iupac.org/home/publications/e-resources/inchi.html
.. _IUPAC: http://www.iupac.org
.. _Tripos: http://www.tripos.com
.. _MOL: http://en.wikipedia.org/wiki/Chemical_table_file
.. _MOL2: http://www.tripos.com/mol2/mol2_format3.html
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

.. _Getting Started:
Getting Started
===============
.. _Installation:

ChemicalToolBoX can be installed on all common unix systems. However, it is developed on Linux and I don't have access to OS X. You are welcome to help improving this documentation, just contact_ me.

.. _contact: https://github.com/bgruening

Prerequisites::

* Python 2.6 or above
* standard C compiler, C++ and Fortran compiler
* Autotools
* CMake
* cairo development files (used for PNG depictions)
* python development files
* Java Runtime Environment (JRE, used by OPSIN and NPLS)

- Debian based systems: apt-get install build-essential gfortran cmake mercurial libcairo2-dev python-dev
- Fedora: yum install make automake gcc gcc-c++ gcc-gfortran cmake mercurial libcairo2-devel python-devel
- OS X (MacPorts_): port install gcc cmake automake mercurial cairo-devel

.. _MacPorts: http://www.macports.org/

1. Clone the latest `Galaxy platform`_::

	hg clone https://https://bitbucket.org/galaxy/galaxy-central/

.. _Galaxy platform: http://wiki.galaxyproject.org/Admin/Get%20Galaxy

2. Navigate to the galaxy-central folder and update it::
	
	cd ~/galaxy-central
	hg pull
	hg update

3. Create folders for toolshed and dependencies::

	mkdir ~/shed_tools
	mkdir ~/galaxy-central/tool_deps

4. Create configuration file::

	cp ~/galaxy-central/universe_wsgi.ini.sample ~/galaxy-central/universe_wsgi.ini

5. Open universe_wsgi.ini and change the dependencies directory::

	LINUX: gedit ~/galaxy-central/universe_wsgi.ini
	OS X: open -a TextEdit ~/galaxy-central/universe_wsgi.ini

6. Search for ``tool_dependency_dir = None`` and change it to ``tool_dependency_dir = ./tool_deps``

7. Remove the hash in front of ``tool_config_file`` and ``tool_path``

8. (Re-)Start the galaxy daemon::

	GALAXY_RUN_ALL=1 sh run.sh --stop-daemon
	GALAXY_RUN_ALL=1 sh run.sh --daemon

After launching galaxy is accessible via the browser at ``http://localhost:8080/``.

To improve the overall performance of NumPy_ you need to disable CPU throttling::

	cpufreq-selector -g performance

.. _NumPy: http://www.numpy.org/



.. _Admin Account
=============
Admin Account
=============

- Register a new account

- Promote user to admin
	- open universe_wsgi.ini
	- search ``admin_users = None`` and change it to ``admin_users = YOUR_EMAIL_ADDRESS``


.. _Toolshed
Toolshed
========

===========================
Installation via webbrowser
===========================

- go to the `admin page`_
- select *Search and browse tool sheds*
- Galaxy test tool shed > Computational chemistry > chemicaltoolbox
- install chemicaltoolbox

.. _admin page: http://localhost:8080/admin


================
API Installation
================

- Generate an `API Key`_
- Run the installation script::
	
	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY 
	-l http://localhost:8080 --url http://testtoolshed.g2.bx.psu.edu/ -o bgruening 
	-r 2c9d1a52824d --name chemicaltoolbox --tool-deps --repository-deps 
	--panel-section-name ChemicalToolBoX

.. _API Key: http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key

========================
JMol Editor Installation
========================

`JMol Editor`_ needs be run on a webserver, this is how to setup the server:

.. _JMol Editor: http://wiki.jmol.org/index.php/Jmol_as_editor

- copy the directory ``jmoleditor`` from /galaxytools/chemicaltoolbox/data_source/ into your Galaxy Root directory ::

	cp -a ~/galaxytools/chemicaltoolbox/data_source/jmoleditor/ ~/galaxy-central/

- launch the webserver ::

	python -m SimpleHTTPServer &

.. _Tools
Tools
=====

- Get Chemical Data
	- JMol Editor
		JMol_ is a viewer of molecular structures but the JMol Editor can be used to alter atom positions or identities and to add/remove atoms.

.. _JMol: http://jmol.sourceforge.net/
	- Online data
		Upload data via FTP or HTTP and load them into your history. Supportes gz/gzip- and rar-files.
	- PubChem download
		Download all molecules from PubChem_ and store them as SMILES file.

.. _PubChem: http://pubchem.ncbi.nlm.nih.gov/

- Chemical Converters
	- Compound converter
		Compound converter joins several `Open Babel command prompt converters`_ in an easy to use tool. It converts various chemistry and moleculare modeling data files. The output format can be specified as well as several parameters. Some parameters are available for all tools (e.g. protonation state & pH) others are specific for a given output format (e.g. exclude isotopes for conversion to canSMI).
	- Molecule recognition
		OSRA_ (Optical Structure Recognition Application) is a utility designed to convert graphical representations of chemical structures into SMILES or SDF. It generates the SMILES or SDF representation of any molecular structure image within a document which is parseable by ImageMagick.
	- IUPAC name-to-structure
		OPSIN_ is a IUPAC name-to-structure conversion tool offering high recall and precision on organic chemical nomenclature.

- Filter / Sort
	- (Multi) Compound search
		Uses the Open Babel Obgrep_ to search for molecules inside multi-molecule files (e.g. SMI, SDF, etc.) or across multiple files.
	- Remove counterions and fragments
		Parses a multiple molecules file and deletes any present counterions or fragments.
	- Remove duplicated molecules
		Filters a library of compounds and removes duplicated molecules comparing either InChI or SMI.
	- Filter
		Filters a library of compounds based on user-defined physico-chemical parameters or predefined options (e.g. Ro5, lead-like properties, etc.). Multiple parameters can be selected for more specific queries. 
	- Remove small molecules
		Filters a library of compounds and removes small molecules below a predefined input number of atoms.

- Search
	- |Spectrophores (TM)| search
		|Spectrophores (TM)| is a screening technology by Silicos_ which converts three-dimensional molecular property data into one-dimensional spectra. Typical characteristics that can be converted include electrostatic potentials, molecular shape, lipophilicity, hardness and softness potentials. The computation is independent of the position and orientation of a molecule and allows an easy comparison of |Spectrophores (TM)| of different molecules.

		Molecules with similar three-dimensional properties and shape, and therefore also similar biological activities, always have similar |Spectrophores (TM)|. As a result this technique is a very powerful tool to investigate the similarity of molecules and can be applied as a screening tool for molecular databases, virtual screening, and database characterisations.
	- Similarity search
		Similarity searches using a variety of different fingerprints using either the chemfp_ FPS type or the `Open Babel` Fastsearch_ index.
	- Substructure search
		Substructure search is based on Open Babel FastSearch_. FastSearch uses molecular fingerprints to prepare and search an index of a multi-molecule datafile.

- Calculate / Modify
	- Compute physico-chemical properties
		Computes several physico-chemical properties (e.g. logP, PSA, MW, etc.) for a set of molecules. Accepts SDF or MOL2 as input file as 3D coordinates of the molecules have to be provided.
	- Add hydrogen atoms
		Parses a molecular file and adds hydrogen atoms at a user-defined pH value.
	- Remove protonation state
		Parses a molecular file and removes the protonation state of every atom.
	- Change title
		Changes the title of a molecule file to a metadata value of a given ID in the same molecule file.
	- Confab
		Confab_ is a conformation generator. The algorithm starts with an input 3D structure which, after some initialisation steps, is used to generate multiple conformers which are filtered on-the-fly to identify diverse low energy conformers.
	- Molecules to fingerprints
		Chemfp_ is a tool for fingerprint generation. It supports the FPS fingerprint file format using `Open Babel`_, OpenEye_ and RDKit_ .
	- SDF to fingerprint
		Read an input SD file (pubchem), extract the fingerprints and store them in a FPS-file.
	- Drug-likeness
		Describes the similarity of a compound to known drugs. Comes with three applicable varieties (QED\ :sub:`w,mo`\ , QED\ :sub:`w,max`\ , QED\ :sub:`w,u` ).
	- Descriptors by RDKit_
		An open source cheminformatics and machine learning toolkit with a lot of overlap with OpenBabel. It therefor can be used to compare results with OpenBabel. The tool offers different descriptor and fingerprint calculations.
	- `Natural Product likeness`_
		Calculates the Natural Product(NP)-likeness of a molecule, i.e. the similarity of the molecule to the structure space covered by known natural products.
	- |Shape-it (TM)|
		|Shape-it (TM)| is a `silicos-it tool`_ that aligns a reference molecule against a set of database molecules using the shape of the molecules as the align criterion. It is based on the use of `gaussian volumes as descriptor for molecular shape`_ as it was introduced by Grant and Pickup.

		|Shape-it (TM)| is a program that is instructed by means of command line options. The program expects a single reference molecule (with three-dimensional coordinates) and a database file containing one or more molecules (with three-dimensional coordinates) that need to be shape-aligned onto the reference molecule. The tool returns all aligned database molecules and their respective shape overlap scores, or the top-best scoring molecules.

	- |Strip-it (TM)|
		|Strip-it (TM)| is a `program by silicos-it`_ that identifies and extracts predefined scaffolds from organic small molecules. The program is linked against the open source C++ library of Open Babel.

		The program comes with a number of predefined molecular scaffolds for extraction. These scaffolds include, amongst others `molecular frameworks`_ as originally described by Bemis and Murcko, `molecular frameworks and the reduced molecular frameworks`_ as described by Ansgar Schuffenhauer and coworkers and `scaffold topologies`_ as described by Sara Pollock and coworkers.

- Chemical Clustering
	- NxN clustering
		Generates hierarchical clusters and visualises clusters with dendrograms. Powered by chemfp_.
	- Taylor-Butina clustering
		`Taylor-Butina clustering`_ is an unsupervised non-hierarchical clustering method which guarantees that every cluster contains molecules which are within a distance cutoff of the central molecule. Powered by chemfp_.

- Fragmentation
	- Fragmenter
		Splits a molecule on predefined spots, e.g. the RECAP-rules.
	- Merging
		Merges small molecules together to larger compounds using  predefined reactions. The options *iteration depth* and *number of repeats* can be used to adjust the created number of compounds and the actual computation time.

- Visualisation
	- Visualisation
		Creates an .svg or .png image of a small set of molecules (few hundreds). Based on `Open Babel`_ PNG_/SVG_ 2D depiction.

.. |Spectrophores (TM)| unicode:: Spectrophores U+2122
.. |Strip-it (TM)| unicode:: Strip-it U+2122
.. |Shape-it (TM)| unicode:: Shape-it U+2122
   .. trademark sign

.. _OPSIN: https://bitbucket.org/dan2097/opsin/overview
.. _program by silicos-it: http://silicos-it.com/software/strip-it/1.0.1/strip-it.html
.. _silicos-it tool: http://silicos-it.com/software/shape-it/1.0.1/shape-it.html
.. _molecular frameworks: http://www.ncbi.nlm.nih.gov/pubmed/8709122
.. _molecular frameworks and the reduced molecular frameworks: http://peter-ertl.com/reprints/Schuffenhauer-JCIM-47-47-2007.pdf
.. _scaffold topologies: http://www.ncbi.nlm.nih.gov/pubmed/18605680
.. _gaussian volumes as descriptor for molecular shape: http://pubs.acs.org/doi/abs/10.1021/j100011a016
.. _obgrep: http://openbabel.org/wiki/Obgrep
.. _FastSearch: http://openbabel.org/wiki/FastSearch
.. _Silicos: http://www.silicos.be/technologies/spectrophore
.. _chemfp: http://chemfp.com/
.. _Open Babel command prompt converters: http://openbabel.org/docs/2.3.0/FileFormats/Overview.html
.. _Open Babel: http://openbabel.org/wiki/Main_Page
.. _OpenEye: http://www.eyesopen.com/
.. _RDKit: http://www.rdkit.org/
.. _Taylor-Butina clustering: http://www.redbrick.dcu.ie/~noel/R_clustering.html
.. _PNG: http://openbabel.org/docs/dev/FileFormats/PNG_2D_depiction.html
.. _SVG: http://openbabel.org/docs/dev/FileFormats/SVG_2D_depiction.html
.. _OSRA: http://cactus.nci.nih.gov/osra/
.. _Confab: https://code.google.com/p/confab/
.. _`natural Product likeness`: http://sourceforge.net/projects/np-likeness/

.. _Bug Tracker
Bug Tracker
===========
Have a bug or a feature request? `Please write a new card`_. Before writing a new card, please search for existing issues.

.. _Please write a new card: https://trello.com/b/t9Wr8lSY

.. _Contributing
Contributing
============
We encourage you to contribute to ChemicalToolBoX! Check out our `Trello board`_ or contact us via e-mail_.

.. _Trello board: https://trello.com/b/t9Wr8lSY
.. _e-mail: bjoern_dot_gruening@gmail.com
