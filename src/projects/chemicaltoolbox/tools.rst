---
title: List of Tools
layout: default
---

=====
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
		Simsearch is a chemfp_ tool searching a FPS file for similar fingerprints.
	- Substructure search
		Substructure search is based on Open Babel FastSearch_. FastSearch uses molecular fingerprints to prepare and search an index of a multi-molecule datafile. It allows fast substructure and structural similarity searching. The indexing is a slow process (~30 minutes for a 250,000 molecule file). The subsequent seaching is much faster, a few seconds, and so can be done interactively.

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
		Read an input SDF file, extract the fingerprints and store them in a fps-file.
	- Drug-likeness
		Describes the similarity of a compound to known drugs. Comes with three applicable varieties (QED\ :sub:`w,mo`\ , QED\ :sub:`w,max`\ , QED\ :sub:`w,u` ).
	- Descriptors by RDKit_
		An open source cheminformatics and machine learning toolkit with a lot of overlap with OpenBabel. It therefor can be used to compare results with OpenBabel. The tool offers different descriptor and fingerprint calculations.
	- Natural Product
		Calculates the Natural Product(NP)-likeness of a molecule, i.e. the similarity of the molecule to the structure space covered by known natural products.
	- |Shape-it (TM)|
		|Shape-it (TM)| is a `silicos-it tool`_ that aligns a reference molecule against a set of database molecules using the shape of the molecules as the align criterion. It is based on the use of `gaussian volumes as descriptor for molecular shape`_ as it was introduced by Grant and Pickup.

		|Shape-it (TM)| is a program that is instructed by means of command line options. The program expects a single reference molecule (with three-dimensional coordinates) and a database file containing one or more molecules (with three-dimensional coordinates) that need to be shape-aligned onto the reference molecule. The tool returns all aligned database molecules and their respective shape overlap scores, or the top-best scoring molecules.

	- |Strip-it (TM)|
		|Strip-it (TM)| is a `program by silicos-it`_ that identifies and extracts predefined scaffolds from organic small molecules. The program is linked against the open source C++ library of Open Babel.

		The program comes with a number of predefined molecular scaffolds for extraction. These scaffolds include, amongst others `molecular frameworks`_ as originally described by Bemis and Murcko, `molecular frameworks and the reduced molecular frameworks`_ as described by Ansgar Schuffenhauer and coworkers and `scaffold topologies`_ as described by Sara Pollock and coworkers.

- Chemical Clustering
	- NxN clustering
		Generates hierarchical clusters and visualises clusters with dendrograms. Accepts fingerprints in FPS format as input.
	- Taylor-Butina clustering
		`Taylor-Butina clustering`_ is an unsupervised non-hierarchical clustering method which guarantees that every cluster contains molecules which are within a distance cutoff of the central molecule. Fingerprints in FPS format are needed as input.

- Fragmentation
	- Fragmenter
		Splits a molecule on predefined spots, following the RECAP-rules.
	- Merging
		Merges small molecules together to larger compounds using  predefined reactions. The options *Molecule dependend iteration depth* and *Number of repeats* can be used to adjust the created number of compounds and the actual computation time.

- Visualisation
	- Visualisation
		Creates an .svg or .png image of a small set of molecules (few hundreds). Based on Open Babel PNG_/SVG_ 2D depiction.

.. |Spectrophores (TM)| unicode:: Spectrophores
.. |Strip-it (TM)| unicode:: Strip-it
.. |Shape-it (TM)| unicode:: Shape-it
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
.. _Chemfp: https://chemfp.readthedocs.org/en/latest/
.. _Open Babel command prompt converters: http://openbabel.org/docs/2.3.0/FileFormats/Overview.html
.. _Open Babel: http://openbabel.org/wiki/Main_Page
.. _OpenEye: http://www.eyesopen.com/
.. _RDKit: http://www.rdkit.org/
.. _Taylor-Butina clustering: http://www.redbrick.dcu.ie/~noel/R_clustering.html
.. _PNG: http://openbabel.org/docs/dev/FileFormats/PNG_2D_depiction.html
.. _SVG: http://openbabel.org/docs/dev/FileFormats/SVG_2D_depiction.html
.. _OSRA: http://cactus.nci.nih.gov/osra/
.. _Confab: https://code.google.com/p/confab/



