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

- `ChemcialToolBoX Installation`_
- Troubleshooting_
- Tools_
- Workflows_
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

CTB is available as a public test instance @ http://ctb.pharmaceutical-bioinformatics.org.


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
	Structure-data-file that can consist of many molecules. Molecules are separated by four Dollar signs ($$$$). Allows the storing of metainformation like molecular mass or a unique identifier. Developed by MDL Information System (Accelrys_).
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

All filetypes are interchangable due to three easy to use converting options:


- the built-in conversion via the pencil icon |pencilicon|

	.. |pencilicon| image:: https://github.com/bgruening/galaxytools/raw/master/chemicaltoolbox/pencil_icon.png
	.. image:: https://github.com/bgruening/galaxytools/raw/master/chemicaltoolbox/convert_pencil.jpg

- the Compound Converter tool described in the Tools section

- the automatic conversion each tool offers

	.. image:: https://github.com/bgruening/galaxytools/raw/master/chemicaltoolbox/internal_conversion.png

.. _Getting Started:

Getting Started
===============

The easiest way to get the ChemicalToolBox up and running is to use CTB in our Galaxy Docker flavor from https://github.com/bgruening/docker-galaxy-chemicaltoolbox. Simply run::

	docker run -i -t -p 8080:80 quay.io/bgruening/galaxy-chemicaltoolbox

And open your webbrowser on `localhost:8080`. You can also install CTB manually, in the following we will guide you through this process.



.. _Installation:

ChemicalToolBoX can be installed on all common Unix systems. 
However, it is developed on Linux and I don't have access to OS X. You are welcome to help improving this documentation, just contact_ me.

For any additional information, especially cluster configuration or general Galaxy_ questions, 
please have a look at the Galaxy Wiki.

- http://wiki.galaxyproject.org/

- http://wiki.galaxyproject.org/Admin/

- http://galaxyproject.org/search/web/

.. _contact: https://github.com/bgruening
.. _Galaxy: http://galaxyproject.org/

Prerequisites::

* Python 2.6 or 2.7
* standard C compiler, C++ and Fortran compiler
* Autotools
* CMake
* cairo development files (used for PNG depictions)
* python development files
* libblas and liblapack development files
* Java Runtime Environment (JRE, used by OPSIN and NPLS)

To install all of the prerequisites you can run the following command, depending on your OS:

- Debian based systems: apt-get install build-essential gfortran cmake mercurial libcairo2-dev python-dev
- Fedora: yum install make automake gcc gcc-c++ gcc-gfortran cmake mercurial libcairo2-devel python-devel
- OS X (MacPorts_): port install gcc cmake automake mercurial cairo-devel

.. _MacPorts: http://www.macports.org/


===================
Galaxy installation
===================


0. Create a sand-boxed Python using virtualenv_ (not necessary but recommended)::

        wget https://raw.github.com/pypa/virtualenv/master/virtualenv.py
	python ./virtualenv.py --no-site-packages galaxy_env
	. ./galaxy_env/bin/activate

.. _virtualenv: http://www.virtualenv.org/


1. Clone the latest `Galaxy platform`_::

	git clone https://github.com/galaxyproject/galaxy.git

.. _Galaxy platform: http://wiki.galaxyproject.org/Admin/Get%20Galaxy

2. Navigate to the galaxy folder and update it::
	
	cd ~/galaxy
	git pull

   This step is not necessary if you have a fresh checkout. Anyway, it is good to know ;)

3. Create folders for toolshed and dependencies::

	mkdir ~/shed_tools
	mkdir ~/galaxy/tool_deps

4. Create configuration file::

	cp ~/galaxy/config/galaxy.ini.sample ~/galaxy/config/galaxy.ini

5. Open config/galaxy.ini and change the dependencies directory::

	LINUX: gedit ~/galaxy/config/galaxy.ini
	OS X: open -a TextEdit ~/galaxy/config/galaxy.ini

6. Search for ``tool_dependency_dir = None`` and change it to ``tool_dependency_dir = ./tool_deps``, remove the ``#`` if needed

7. Remove the ``#`` in front of ``tool_config_file`` and ``tool_path``

8. (Re-)Start the galaxy daemon::

	sh run.sh --reload
	
   In deamon mode all logs will be written to main.log in your Galaxy Home directory. You can also use::
   
	run.sh   

   During the first startup Galaxy will prepare your database. That can take some time. Have a look at the log file if you want to know what happens.

After launching galaxy is accessible via the browser at ``http://localhost:8080/``.


.. _Admin Account:

=======================
Tool Shed configuration
=======================

- Register a new user account in your Galaxy instance: Top Panel → User → Register
- Become an admin
	- open ``config/galaxy.ini`` in your favourite text editor (gedit config/galaxy.ini)
	- search ``admin_users = None`` and change it to ``admin_users = EMAIL_ADDRESS`` (your Galaxy Username)
	- remove the ``#`` if needed
- restart Galaxy

::

	sh run.sh --reload

.. _ChemcialtoolboX Installation:

============================
ChemicalToolBoX installation
============================

ChemicalToolBoX will automatically download and compile all requirements, 
like `Open Babel`_, RDKit_, chemfp_, numpy_ and so on. It can take up to 2-3 hours.


Installation via Galaxy API (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Generate an `API Key`_
- Run the installation script::
	
	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o bgruening --name chemicaltoolbox --tool-deps --repository-deps --panel-section-name ChemicalToolBoX

You can watch the installation status under: Top Panel → Admin → Manage installed tool shed repositories


.. _API Key: http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key
.. _`test tool shed`: http://testtoolshed.g2.bx.psu.edu/


Installation via webbrowser
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- go to the `admin page`_
- select *Search and browse tool sheds*
- Galaxy tool shed > Computational chemistry > chemicaltoolbox
- install chemicaltoolbox

.. _admin page: http://localhost:8080/admin


Additional Notes
~~~~~~~~~~~~~~~~

You can also configure CTB to use system installed binaries, but you will loose some degree of reproducibility.
Nevertheless, if you want to do this the recommended depencency versions are specified in a file called
``tool_dependencies.xml``, located in each subfolder.





.. _Troubleshooting:

===============
Troubleshooting
===============

If you have any trouble or the installation did not finish properly, do not hesitate to contact me. However, if the 
installation fails during the Galaxy installation, you can have a look at the `Galaxy wiki`_. If the ChemicalToolBoX installation fails, 
you can try to run::

	python ./scripts/api/repair_tool_shed_repository.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o bgruening -r 30ae0e5218b4 --name chemicaltoolbox

That will rerun all failed installation routines. Alternatively, you can navigate to the ChemicalToolBoX repository in 
your browser and repair manually: 
Top Panel → Admin → Manage installed tool shed repositories → chemicaltoolbox → Repository Actions → Repair repository

------


On slow computers and during the compilation of large software libraries, like openbabel or boost, 
the Tool Shed can run into a timeout and kills the installation.
That problem is known and should be fixed in the near future.

If you encouter a timeout or 'hung' during the installation you can increase the ``threadpool_kill_thread_limit`` in your `config/galaxy.ini` file.


------

**Database locking errors**

Please note that Galaxy per default uses a SQLite database. Sqlite is not intended for production use. 
With multiple users or complex components, like that workflow, you will see database locking errors. 
We highly recommend to use PostgreSQL for any kind of production system.


.. _Galaxy wiki: http://wiki.galaxyproject.org/


========================
Jmol Editor Installation
========================

`Jmol Editor`_ needs be run on a separate webserver, this is how to setup the server:

.. _Jmol Editor: http://wiki.jmol.org/index.php/Jmol_as_editor


- download Jmol Editor from::

	wget https://github.com/bgruening/download_store/raw/master/jmoleditor.tar.gz

- copy the directory ``jmoleditor`` into your Galaxy Root directory ::

	cp -a ~/galaxytools/chemicaltoolbox/data_source/jmoleditor/ ~/galaxy/

- launch the webserver from your galaxy root directory ::

	python -m SimpleHTTPServer &

.. _Tools:

Tools
=====

- Get Chemical Data
	- Jmol Editor
		Jmol_ Editor can be used to paint structures or alter atoms or identities from single molecules.

.. _Jmol: http://jmol.sourceforge.net/
	- Online data
		Upload data via FTP or HTTP and load them into your history. Supportes compressed-files.
	- PubChem download
		Download all molecules from PubChem_ and store them in a single large SMILES file.

.. _PubChem: http://pubchem.ncbi.nlm.nih.gov/

- Chemical Converters
	- Compound converter
		Compound converter joins several `Open Babel command prompt converters`_ in an easy to use tool. It converts various chemistry and moleculare modeling data files. The output format can be specified as well as several parameters. Some parameters are available for all tools (e.g. protonation state & pH) others are specific for a given output format (e.g. exclude isotopes for conversion to canonical SMILES).
	- Molecule recognition
		OSRA_ (Optical Structure Recognition Application) is a utility designed to convert graphical representations of chemical structures into SMILES or SDF. It generates the SMILES or SDF representation of any molecular structure image within a document which is parseable by GraphicMagick.
	- IUPAC name-to-structure
		OPSIN_ is a IUPAC name-to-structure conversion tool offering high recall and precision on organic chemical nomenclature.

- Filter / Sort
	- (Multi) Compound search
		Uses the Open Babel Obgrep_ to search for molecules inside multi-molecule files (e.g. SMI, SDF, etc.).
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
		10 different fingerprints can be calculated from all common file formats using chemfp_. Chemfp supports the FPS fingerprint file format and is utilising `Open Babel`_, OpenEye_ and RDKit_.
	- SDF to fingerprint
		Read an input SD file (PubChem), extract the fingerprints and store them in a FPS-file.
	- Drug-likeness
		Estimates the drug-likeness of molecules and reports a score. Comes with three applicable varieties (QED\ :sub:`w,mo`\ , QED\ :sub:`w,max`\ , QED\ :sub:`w,u` ).
	- Descriptors by RDKit_
		This tool calculates all available descriptors from RDKit_..
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
	- Depiction
		Creates an .svg or .png image of a small set of molecules (few hundreds). Based on `Open Babel`_ PNG_/SVG_ 2D depiction.
	- More to come ...
		We are working on several ideas how to improve the visualision of small and large libraries in Galaxy. If
		you are interested and want to discuss it further please contact me (e-mail_).


.. _Workflows:

Workflows
=========

An example workflow is located in the `Tool Shed`::

	 http://toolshed.g2.bx.psu.edu/view/bgruening/chemicaltoolbox_merging_chemical_databases_workflow 

You can install the workflow with the API::

	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o bgruening -r e1bc8415f875 --name chemicaltoolbox_merging_chemical_databases_workflow --tool-deps --repository-deps --panel-section-name ChemicalToolBoX

or as described above via webbrowser. You have now successfully installed the workflow, 
to import it to all your users you need to go to the admin panel, choose the worklow and import it.
For more information have a look at the Galaxy wiki::

	http://wiki.galaxyproject.org/ToolShedWorkflowSharing#Finding_workflows_in_tool_shed_repositories

Please **note** that Galaxy per default uses a SQLite database. Sqlite is not intended for production use. 
With multiple users or complex components, like that workflow, you will see database locking errors. 
We highly recommend to use PostgreSQL for any kind of production system.



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
.. _numpy: http://www.numpy.org/
.. _`natural Product likeness`: http://sourceforge.net/projects/np-likeness/

.. _Bug Tracker:


Publications using CTB
======================

 * `The Purchasable Chemical Space: a Detailed Picture <http://pubs.acs.org/doi/abs/10.1021/acs.jcim.5b00116>`_


Bug Tracker
===========
Have a bug or a feature request? `Please write a new card`_. Before writing a new card, please search for existing issues.

.. _Please write a new card: https://trello.com/b/t9Wr8lSY/chemicaltoolbox

.. _Contributing:

Contributing
============
We encourage you to contribute to ChemicalToolBoX! Check out our `Trello board`_ or contact us via e-mail_.

.. _Trello board: https://trello.com/b/t9Wr8lSY/chemicaltoolbox
.. _e-mail: bjoern.gruening@gmail.com
