---
title: Installation instruction
layout: default
---

===============
Getting Started
===============

ChemicalToolBoX can be installed on all common unix systems. 
However, it is developed on Linux and I don't have access to OS X. You are welcome to help improving this documentation, just contact_ me.

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


=============
Admin Account
=============

- Register a new account

- Promote user to admin
	- open universe_wsgi.ini
	- search ``admin_users = None`` and change it to ``admin_users = YOUR_EMAIL_ADDRESS``


========
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



