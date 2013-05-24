---
title: Installation instruction
layout: default
---

===============
Getting Started
===============

ChemicalToolBoX can be installed on all common unix systems.
However, it is developed on Linux and I don't have access to OSX.
If you want to help to enhance the documentation, please contact_ me.

.. _contact: https://github.com/bgruening

Prerequisites::

* Python 2.6 or above
* standard C, C++ and Fortran compiler
* Autotools
* CMake

    - Debian based systems: apt-get install build-essential gfortran cmake mercurial
    - Fedora: yum install make automake gcc gcc-c++ gcc-gfortran cmake mercurial
    - OSX (MacPorts_): port install gcc cmake automake mercurial

.. _macports: http://www.macports.org/


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

4. Open universe_wsgi.ini and change the dependencies directory::

	gedit ~/galaxy-central/universe_wsgi.ini


5. Search for ``tool_dependency_dir = None`` and change it to ``tool_dependency_dir = ./tool_deps``

6. (Re-)Start the galaxy daemon::

	GALAXY_RUN_ALL=1 sh run.sh --stop-daemon
	GALAXY_RUN_ALL=1 sh run.sh --daemon

After launching galaxy is accessible via the browser at ``http://localhost:8080/``.

Admin Account
=============

1. Register a new account

2. Promote user to admin
	- open universe_wsgi.ini
	- search ``admin_users = None`` and change it to ``admin_users = YOUR_EMAIL_ADDRESS``

Toolshed
========

**Installation via webbrowser**

- go to the `admin page`_
- select *Search and browse tool sheds*
- Galaxy test tool shed >> Computational chemistry >> chemicaltoolbox
- install chemicaltoolbox

**Installtion via the shell**

- generate an API_ Key
- execute the following code from your terminal

.. code-block:: console

   $ ./scripts/api/install_tools.py --api  <Your Galaxy API Key> \
    -l http://localhost:8080 --url http://testtoolshed.g2.bx.psu.edu/ \
    -o bgruening -r 7e98219aa915  --name chemicaltoolbox \
    --tool-deps --repository-deps \
    --panel-section-name ChemicalToolBoX


.. _admin page: http://localhost:8080/admin
.. _API: http://wiki.galaxyproject.org/Admin/API#Generate_the_Admin_Account_API_Key



