This package is a Galaxy workflow for BlockClust pipeline.


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



.. _Getting Started:

Getting Started
===============

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

	hg clone https://bitbucket.org/galaxy/galaxy-central/

.. _Galaxy platform: http://wiki.galaxyproject.org/Admin/Get%20Galaxy

2. Navigate to the galaxy-central folder and update it::
	
	cd ~/galaxy-central
	hg pull
	hg update
   
   This step is not necessary if you have a fresh checkout. Anyway, it is good to know ;)

3. Create folders for toolshed and dependencies::

	mkdir ~/shed_tools
	mkdir ~/galaxy-central/tool_deps

4. Create configuration file::

	cp ~/galaxy-central/universe_wsgi.ini.sample ~/galaxy-central/universe_wsgi.ini

5. Open universe_wsgi.ini and change the dependencies directory::

	LINUX: gedit ~/galaxy-central/universe_wsgi.ini
	OS X: open -a TextEdit ~/galaxy-central/universe_wsgi.ini

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
	- open ``universe_wsgi.ini`` in your favourite text editor (gedit universe_wsgi.ini)
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
	
	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o rnateam -r e9b2400cc569 --name blockclust_workflow --tool-deps --repository-deps --panel-section-name ChemicalToolBoX

The -r argument specifies the version of ChemicalToolBoX. You can get the latest revsion number from the 
`test tool shed`_ or with the following command::

	hg identify http://toolshed.g2.bx.psu.edu/repos/bgruening/chemicaltoolbox

You can watch the installation status under: Top Panel → Admin → Manage installed tool shed repositories


.. _API Key: http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key
.. _`test tool shed`: http://testtoolshed.g2.bx.psu.edu/


Installation via webbrowser
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- go to the `admin page`_
- select *Search and browse tool sheds*
- Galaxy test tool shed > Sequence Analysis  > blockclust_workflow
- install chemicaltoolbox

.. _admin page: http://localhost:8080/admin



.. _Troubleshooting:

===============
Troubleshooting
===============

If you have any trouble or the installation did not finish properly, do not hesitate to contact me. However, if the 
installation fails during the Galaxy installation, you can have a look at the `Galaxy wiki`_. If the ChemicalToolBoX installation fails, 
you can try to run::

	python ./scripts/api/repair_tool_shed_repository.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o rnateam -r e9b2400cc569 --name blockclust_workflow

That will rerun all failed installation routines. Alternatively, you can navigate to the ChemicalToolBoX repository in 
your browser and repair manually: 
Top Panel → Admin → Manage installed tool shed repositories → chemicaltoolbox → Repository Actions → Repair repository

------


On slow computers and during the compilation of large software libraries, like R, 
the Tool Shed can run into a timeout and kills the installation.
That problem is known and should be fixed in the near future.

If you encouter a timeout or 'hung' during the installation you can increase the ``threadpool_kill_thread_limit`` in your universe_wsgi.ini file.


------

**Database locking errors**

Please note that Galaxy per default uses a SQLite database. Sqlite is not intended for production use. 
With multiple users or complex components, like that workflow, you will see database locking errors. 
We highly recommend to use PostgreSQL for any kind of production system.


.. _Galaxy wiki: http://wiki.galaxyproject.org/


Workflows
=========

An example workflow is located in the `Tool Shed`::

	  http://testtoolshed.g2.bx.psu.edu/view/rnateam/blockclust_workflow

You can install the workflow with the API::

	python ./scripts/api/install_tool_shed_repositories.py --api YOUR_API_KEY -l http://localhost:8080 --url http://toolshed.g2.bx.psu.edu/ -o rnateam -r e9b2400cc569 --name blockclust_workflow --tool-deps --repository-deps --panel-section-name BlockClust

or as described above via webbrowser. You have now successfully installed the workflow, 
to import it to all your users you need to go to the admin panel, choose the worklow and import it.
For more information have a look at the Galaxy wiki::

	http://wiki.galaxyproject.org/ToolShedWorkflowSharing#Finding_workflows_in_tool_shed_repositories

Please **note** that Galaxy per default uses a SQLite database. Sqlite is not intended for production use. 
With multiple users or complex components, like that workflow, you will see database locking errors. 
We highly recommend to use PostgreSQL for any kind of production system.





Sample Data
===========

As an example, we will use the first public assembly of the 2011 Shiga-toxin
producing *Escherichia coli* O104:H4 outbreak in Germany. This was part of the
open-source crowd-sourcing analysis described in Rohde et al. (2011) and here:
https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/wiki

You can upload this assembly directly into Galaxy using the "Upload File" tool
with either of these URLs - Galaxy should recognise this is a FASTA file with
3,057 sequences:

* http://static.xbase.ac.uk/files/results/nick/TY2482/TY2482.fasta.txt
* https://github.com/ehec-outbreak-crowdsourced/BGI-data-analysis/blob/master/strains/TY2482/seqProject/BGI/assemblies/NickLoman/TY2482.fasta.txt

This FASTA file ``TY2482.fasta.txt`` was the initial TY-2482 strain assembled
by Nick Loman from 5 runs of Ion Torrent data released by the BGI, using the
MIRA 3.2 assembler. It was initially released via his blog,
http://pathogenomics.bham.ac.uk/blog/2011/06/ehec-genome-assembly/

We will also need a training set of known *E. coli* genes, for example the
model strain *Escherichia coli* str. K-12 substr. MG1655 which is well
annotated. You can upload the NCBI FASTA file ``NC_000913.ffn`` of the
gene nucleotide sequences directly into Galaxy via this URL, which Galaxy
should recognise as a FASTA file with 4,321 sequences:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.ffn

Then run the workflow, which should produce 2,333 predicted genes for the
TY2482 assembly (two FASTA files, nucleotide and protein sequences).


Citation
========

If you use this workflow directly, or a derivative of it, or the associated
wrappers for Galaxy, in work leading to a scientific publication,
please cite:

P. Videm  at al...

For Glimmer3 please cite:

Delcher, A.L., Bratke, K.A., Powers, E.C., and Salzberg, S.L. (2007)
Identifying bacterial genes and endosymbiont DNA with Glimmer.
Bioinformatics 23(6), 673-679.
http://dx.doi.org/10.1093/bioinformatics/btm009

For EMBOSS please cite:

Rice, P., Longden, I. and Bleasby, A. (2000)
EMBOSS: The European Molecular Biology Open Software Suite
Trends in Genetics 16(6), 276-277.
http://dx.doi.org/10.1016/S0168-9525(00)02024-2


Additional References
=====================

Rohde, H., Qin, J., Cui, Y., Li, D., Loman, N.J., et al. (2011)
Open-source genomic analysis of shiga-toxin-producing E. coli O104:H4.
New England Journal of Medicine 365, 718-724.
http://dx.doi.org/10.1056/NEJMoa1107643


Availability
============

This workflow is available on the main Galaxy Tool Shed:

 http://testtoolshed.g2.bx.psu.edu/view/rnateam/blockclust_workflow 

Development is being done on github:

https://github.com/bgruening/galaxytools/tree/master/workflows/blockclust


Dependencies
============

These dependencies should be resolved automatically via the Galaxy Tool Shed:

* http://testtoolshed.g2.bx.psu.edu/view/iuc/package_samtools_0_1_19 
* http://testtoolshed.g2.bx.psu.edu/view/iuc/package_r_3_0_1
* http://testtoolshed.g2.bx.psu.edu/view/rnateam/package_segemehl_0_1_6 
* http://testtoolshed.g2.bx.psu.edu/view/iuc/msa_datatypes 
* http://testtoolshed.g2.bx.psu.edu/view/iuc/package_infernal_1_1rc4 
* http://testtoolshed.g2.bx.psu.edu/view/rnateam/blockbuster 
* http://testtoolshed.g2.bx.psu.edu/view/bgruening/package_eden_1_1
* http://testtoolshed.g2.bx.psu.edu/view/iuc/package_mcl_12_135 
